// Standalone reproducer for the libstdc++-specific RunStorage read/write-mixing failure that shows up
// as "RunStorage::get_*() stream not good" during the mou chance/stack processing.  No PANTHER, no
// timing: it just drives RunStorage's exact read<->write interleaving sequentially, across many
// orderings.  Build + run on Linux (libstdc++) where it should reproduce deterministically; on mac
// (libc++) it passes.  Usage: rs_rwtest [niter]
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include "RunStorage.h"
using namespace std;

static vector<char> make_serial(int npar, int nobs, double base)
{
    vector<char> v((npar + nobs) * sizeof(double), 0);
    double *d = reinterpret_cast<double *>(v.data());
    for (int i = 0; i < npar + nobs; i++) d[i] = base + i;
    return v;
}

int main(int argc, char **argv)
{
    const int npop = 10, nstack = 100, total = npop + nstack; // mimic chance: 10 pop + 100 nested stack
    const int npar = 6, nobs = 4;
    const int niter = (argc > 1) ? atoi(argv[1]) : 500;
    vector<string> par_names, obs_names;
    for (int i = 0; i < npar; i++) par_names.push_back("par" + to_string(i));
    for (int i = 0; i < nobs; i++) obs_names.push_back("obs" + to_string(i));

    for (int it = 0; it < niter; it++)
    {
        try
        {
            RunStorage rs("rs_rwtest.rns");
            rs.reset(par_names, obs_names, "rs_rwtest.rns");
            vector<int> rids;
            vector<double> mp(npar, 1.0);
            for (int i = 0; i < total; i++)
                rids.push_back(rs.add_run(mp, "GEN=0_MEMBER=" + to_string(i % npop), 0.0));

            // RUN PHASE: interleave dispatch-reads (get_serial_pars/get_info) with result-writes
            // (update_run) in a per-iteration-varied order (stands in for agent-completion order).
            for (int i = 0; i < total; i++)
            {
                int r = rids[(i * 37 + it * 13) % total];
                vector<char> sp = rs.get_serial_pars(r);                 // read (dispatch)
                int st; string itx; double iv;
                rs.get_info(r, st, itx, iv);                             // read (dispatch)
                rs.update_run(r, make_serial(npar, nobs, r * 100.0));    // write (store result)
            }

            // PROCESSING / read-back: mirrors get_run_info_map + get_failed_run_ids + update_from_runs
            for (int i = 0; i < total; i++) { int st; string itx; double iv; rs.get_info(i, st, itx, iv); }
            for (int i = 0; i < total; i++) rs.get_run_status(i);
            for (int i = 0; i < total; i++)
            {
                vector<double> pv, ov; string itx; double iv;
                rs.get_run(i, pv, ov, itx, iv);
            }
            rs.free_memory();
        }
        catch (const exception &e)
        {
            cout << "REPRODUCED at iter " << it << ":  " << e.what() << endl;
            return 1;
        }
    }
    cout << "no failure in " << niter << " iters (this build's stdlib tolerates the read/write mixing)" << endl;
    return 0;
}
