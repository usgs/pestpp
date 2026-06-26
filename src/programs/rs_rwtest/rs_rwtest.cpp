// Standalone reproducer for the libstdc++-specific RunStorage read/write-mixing failure that shows up
// as "RunStorage::get_*() stream not good" during the mou chance/stack processing.  No PANTHER, no
// timing: it drives RunStorage's exact read<->write interleaving sequentially, with a RANDOM dispatch/
// completion order each iteration (mimicking nondeterministic agent ordering).  Build + run on Linux
// (libstdc++) where it should reproduce; on mac (libc++) it passes.
//   Usage: rs_rwtest [niter]            run iters 0..niter-1 (each seeded by its index)
//          rs_rwtest 1 <seed>           replay a single ordering (e.g. the one that reproduced)
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "RunStorage.h"
using namespace std;

static vector<char> make_serial(int npar, int nobs, double base)
{
    vector<char> v((npar + nobs) * sizeof(double), 0);
    double *d = reinterpret_cast<double *>(v.data());
    for (int i = 0; i < npar + nobs; i++) d[i] = base + i;
    return v;
}

static bool run_once(unsigned seed, int total, int npop, int npar, int nobs,
                     const vector<string> &par_names, const vector<string> &obs_names)
{
    std::mt19937 rng(seed);
    RunStorage rs("rs_rwtest.rns");
    rs.reset(par_names, obs_names, "rs_rwtest.rns");
    vector<int> rids;
    vector<double> mp(npar, 1.0);
    for (int i = 0; i < total; i++)
        rids.push_back(rs.add_run(mp, "GEN=0_MEMBER=" + to_string(i % npop), 0.0));

    // RUN PHASE: randomly interleave dispatch-reads (get_serial_pars/get_info) and result-writes
    // (update_run).  at each step randomly either dispatch the next queued run or store a random
    // already-dispatched one - so reads and writes alternate in an unpredictable order.
    vector<int> to_dispatch(rids);
    std::shuffle(to_dispatch.begin(), to_dispatch.end(), rng);
    vector<int> dispatched;
    size_t di = 0;
    while (di < to_dispatch.size() || !dispatched.empty())
    {
        bool can_dispatch = di < to_dispatch.size();
        bool do_store = !dispatched.empty() && (!can_dispatch || (rng() & 1u));
        if (do_store)
        {
            std::uniform_int_distribution<size_t> pick(0, dispatched.size() - 1);
            size_t k = pick(rng);
            int r = dispatched[k];
            dispatched[k] = dispatched.back();
            dispatched.pop_back();
            rs.update_run(r, make_serial(npar, nobs, r * 100.0)); // write
        }
        else
        {
            int r = to_dispatch[di++];
            vector<char> sp = rs.get_serial_pars(r);              // read
            int st; string itx; double iv;
            rs.get_info(r, st, itx, iv);                          // read
            dispatched.push_back(r);
        }
    }

    // PROCESSING / read-back in random order (mirrors get_run_info_map + get_failed_run_ids +
    // update_from_runs over the completed runs).
    vector<int> order(total);
    for (int i = 0; i < total; i++) order[i] = i;
    std::shuffle(order.begin(), order.end(), rng);
    for (int i : order) { int st; string itx; double iv; rs.get_info(i, st, itx, iv); }
    std::shuffle(order.begin(), order.end(), rng);
    for (int i : order) rs.get_run_status(i);
    std::shuffle(order.begin(), order.end(), rng);
    for (int i : order)
    {
        vector<double> pv, ov; string itx; double iv;
        rs.get_run(i, pv, ov, itx, iv);
    }
    rs.free_memory();
    return true;
}

int main(int argc, char **argv)
{
    const int npop = 10, nstack = 100, total = npop + nstack; // mimic chance: 10 pop + 100 nested stack
    const int npar = 6, nobs = 4;
    const int niter = (argc > 1) ? atoi(argv[1]) : 500;
    vector<string> par_names, obs_names;
    for (int i = 0; i < npar; i++) par_names.push_back("par" + to_string(i));
    for (int i = 0; i < nobs; i++) obs_names.push_back("obs" + to_string(i));

    // single-seed replay mode: rs_rwtest 1 <seed>
    if (argc > 2)
    {
        unsigned seed = (unsigned)strtoul(argv[2], nullptr, 10);
        try { run_once(seed, total, npop, npar, nobs, par_names, obs_names); }
        catch (const exception &e) { cout << "REPRODUCED seed " << seed << ":  " << e.what() << endl; return 1; }
        cout << "seed " << seed << " ok" << endl;
        return 0;
    }

    for (int it = 0; it < niter; it++)
    {
        try { run_once((unsigned)it, total, npop, npar, nobs, par_names, obs_names); }
        catch (const exception &e)
        {
            cout << "REPRODUCED at iter/seed " << it << ":  " << e.what()
                 << "   (replay with: rs_rwtest 1 " << it << ")" << endl;
            return 1;
        }
    }
    cout << "no failure in " << niter << " iters (this build's stdlib tolerates the read/write mixing)" << endl;
    return 0;
}
