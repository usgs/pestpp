#!/usr/bin/env python3
"""Guard random-number portability: no implementation-defined RNG in the shared algorithm code.

std::mt19937 is specified bit-for-bit by the standard - the same seed gives the same integer
stream on every platform and every standard library. NOTHING ELSE IN <random> IS. The
distributions and std::shuffle are implementation-defined, so libc++ (macOS), libstdc++ (Linux)
and the MSVC CRT are each free to turn identical engine bits into different values, and they do.

That is not a last-bit difference. std::shuffle on identical engine state yields a COMPLETELY
different permutation, which is how pestpp-mou came to be irreproducible across platforms while
pestpp-ies - whose subset draw goes through the in-house uniform_int_draws() - was fine.

Use instead, all of which are plain arithmetic over the raw engine stream and therefore portable:

    uniform_draws()       uniform_int_draws()      portable_shuffle()

The rule is easy to state and impossible to remember, which is why it is checked rather than
documented. Run from the repo root; wired into CI beside check_public_surface.py.
"""
import os
import re
import sys

# The algorithm code whose results are compared across platforms. This is now the WHOLE of
# pestpp_common: the one file that used to be exempt, DifferentialEvolution.cpp, held the
# deprecated global-optimization DE and has been removed - pestpp-mou's 'de' generator is the
# maintained implementation, and it goes through portable_shuffle().
SCAN_DIRS = ["src/libs/pestpp_common"]
SKIP_FILES = {
    "RedSVD-h.h",          # already portable by construction; carries its own explanation
}

BANNED = [
    (re.compile(r"\bstd::shuffle\s*\("), "std::shuffle"),
    # bare shuffle( too - 'using namespace std' is pervasive here, so the std:: is often absent
    (re.compile(r"(?<![\w:.])shuffle\s*\("), "shuffle"),
    (re.compile(r"\b(?:std::)?uniform_int_distribution\b"), "uniform_int_distribution"),
    (re.compile(r"\b(?:std::)?uniform_real_distribution\b"), "uniform_real_distribution"),
    (re.compile(r"\b(?:std::)?normal_distribution\b"), "normal_distribution"),
    (re.compile(r"\b(?:std::)?random_device\b"), "random_device"),
]

ALLOW = re.compile(r"portable_shuffle|random_shuffle_is_banned")


def strip_comment(line):
    """crude but sufficient: drop // comments so prose about the rule does not trip it."""
    i = line.find("//")
    return line if i < 0 else line[:i]


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    problems = []
    for d in SCAN_DIRS:
        full = os.path.join(root, d)
        for dirpath, _, fnames in os.walk(full):
            for fn in sorted(fnames):
                if not fn.endswith((".cpp", ".h")) or fn in SKIP_FILES:
                    continue
                path = os.path.join(dirpath, fn)
                with open(path, encoding="utf-8", errors="replace") as f:
                    in_block_comment = False
                    for n, raw in enumerate(f, 1):
                        line = raw
                        # skip block comments, where the rule itself is explained
                        if in_block_comment:
                            if "*/" in line:
                                in_block_comment = False
                            continue
                        if "/*" in line and "*/" not in line:
                            in_block_comment = True
                            continue
                        line = strip_comment(line)
                        if ALLOW.search(line):
                            continue
                        for rx, name in BANNED:
                            if rx.search(line):
                                problems.append((os.path.relpath(path, root), n, name,
                                                 line.strip()[:90]))
                                break
    if problems:
        print("RNG PORTABILITY CHECK FAILED\n")
        print("These are implementation-defined: the same seed produces DIFFERENT values on")
        print("macOS, Linux and Windows, so any result derived from them is not reproducible")
        print("across platforms.  Use uniform_draws(), uniform_int_draws() or")
        print("portable_shuffle() - arithmetic over the raw mt19937 stream - instead.\n")
        for path, n, name, text in problems:
            print("  {0}:{1}: {2}".format(path, n, name))
            print("      {0}".format(text))
        print("\n{0} occurrence(s)".format(len(problems)))
        return 1
    print("rng portability check: PASS (no implementation-defined RNG in {0})".format(
        ", ".join(SCAN_DIRS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
