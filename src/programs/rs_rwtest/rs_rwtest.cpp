// Standalone reproducer for the libstdc++-specific RunStorage read/write-mixing failure ("RunStorage::
// get_*() stream not good") seen during mou chance/stack processing.  Sets up one storage of 110 runs,
// then hammers a tight loop of RANDOM read/write operations on it - aggressively exercising the
// write<->read transitions on the single shared fstream, with NO per-step setup cost.  Build + run on
// Linux (libstdc++) where it should reproduce; on mac (libc++) it passes.
//   Usage: rs_rwtest [nops] [seed]      (default 5,000,000 ops, seed 1)
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <random>
#include "RunStorage.h"
using namespace std;

int main(int argc, char **argv)
{
    const int npop = 10, nstack = 100, total = npop + nstack;
    const int npar = 6, nobs = 4;
    const long nops = (argc > 1) ? atol(argv[1]) : 5000000L;
    const unsigned seed = (argc > 2) ? (unsigned)strtoul(argv[2], nullptr, 10) : 1u;

    vector<string> par_names, obs_names;
    for (int i = 0; i < npar; i++) par_names.push_back("par" + to_string(i));
    for (int i = 0; i < nobs; i++) obs_names.push_back("obs" + to_string(i));

    RunStorage rs("rs_rwtest.rns");
    rs.reset(par_names, obs_names, "rs_rwtest.rns");
    vector<double> mp(npar, 1.0);
    for (int i = 0; i < total; i++) rs.add_run(mp, "GEN=0_MEMBER=" + to_string(i % npop), 0.0);

    // a serialized par+obs result buffer reused for update_run writes
    vector<char> serial((npar + nobs) * sizeof(double), 0);

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> pick_run(0, total - 1);
    std::uniform_int_distribution<int> pick_op(0, 5);  // 1 write op, 5 read ops

    long writes = 0, reads = 0;
    for (long n = 0; n < nops; n++)
    {
        if (n % 500000 == 0)
            cerr << "\r  " << n / 1000 << "k / " << nops / 1000 << "k ops  (w=" << writes
                 << " r=" << reads << ")   " << flush;
        int r = pick_run(rng);
        int op = pick_op(rng);
        try
        {
            switch (op)
            {
            case 0: rs.update_run(r, serial); writes++; break;                                  // write
            default:
            {
                reads++;
                int st; string itx; double iv;
                if (op == 1) rs.get_info(r, st, itx, iv);
                else if (op == 2) { vector<double> pv, ov; rs.get_run(r, pv, ov, itx, iv); }
                else if (op == 3) rs.get_run_status(r);
                else if (op == 4) rs.get_serial_pars(r);
                else (void)rs.get_nruns();
                break;
            }
            }
        }
        catch (const exception &e)
        {
            cout << "REPRODUCED after " << n << " ops (writes=" << writes << " reads=" << reads
                 << ", seed=" << seed << "):  " << e.what() << endl;
            return 1;
        }
    }
    cout << "no failure in " << nops << " ops (writes=" << writes << " reads=" << reads
         << "): this build's stdlib tolerates the read/write mixing" << endl;
    rs.free_memory();
    return 0;
}
