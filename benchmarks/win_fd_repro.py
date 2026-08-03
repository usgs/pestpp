"""Hunt the windows file-descriptor corruption.

Twice now, a run through the C API on windows has ended with a MODEL's console output written
INSIDE one of pest++'s own output files - once in pest.phi.meas.csv, once in
pest.phi.composite.csv. Both times the only symptom was "1 of 25 files differ", which reads
like a numeric wobble; it is not, and the second occurrence is what proved it.

Why those files and not others: the FileManager holds the phi csvs (and the .rec) open for the
whole run, and RunStorage holds its .rns fstream open too, so those descriptors are live across
every model spawn. The ensemble csvs are opened, written and closed between runs and have never
been seen corrupted.

It is intermittent, so this loops. Run it on windows from the benchmarks directory:

    python win_fd_repro.py            # 30 iterations
    python win_fd_repro.py 100        # more

A HIT prints the offending file and the first lines of model output found inside it. Send that,
plus the master's .rec and .rmr from the named directory, and the run is diagnosable.

Clean runs are weak evidence - the race may simply not have fired - so prefer a long loop.
"""
import os
import sys
import traceback

import api_vs_exe_tests as t


def main(n_iter):
    hits = 0
    for i in range(n_iter):
        try:
            t.api_vs_exe_ies_test()
        except AssertionError as e:
            msg = str(e)
            hits += 1
            print("\n=== iteration {0}: FAILURE ===".format(i))
            print(msg[:2000])
            # the detector names corruption explicitly; anything else is a real difference
            if "corrupted" in msg or "model console output" in msg:
                print("\n>>> THIS IS THE DESCRIPTOR CORRUPTION <<<")
            else:
                print("\n>>> a VALUE difference, not corruption - also worth reporting <<<")
            break
        except Exception:
            print("\n=== iteration {0}: unexpected error ===".format(i))
            traceback.print_exc()
            break
        else:
            print("iteration {0}: clean".format(i), flush=True)

    # whether or not the comparison failed, scan both directories directly: corruption could
    # land in a file the comparison does not flag, which would otherwise pass silently
    for tag in ("apiexe_ies_exe", "apiexe_ies_api"):
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), tag)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if not f.endswith(".csv"):
                continue
            why = t.looks_corrupted(os.path.join(d, f))
            if why:
                hits += 1
                print("\n>>> CORRUPTION in {0}/{1}: {2}".format(tag, f, why))

    print("\n{0} iteration(s) run, {1} hit(s)".format(n_iter, hits))
    return 1 if hits else 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    sys.exit(main(n))
