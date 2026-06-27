// Targets the MOU *nested chance* queue->run->process pattern against RunStorage, to chase the
// chance-only "RunStorage::get_*() stream not good" failure (generic random read/write access passes,
// so the bug is specific to this pattern).  Faithfully mimics MOEA::run_population +
// Constraints::add_runs + process_stack_runs/update_from_runs:
//   1. reset the storage (run manager reinitialize)
//   2. QUEUE: stack runs first, per dv member (run_ids 0..nstack-1), recording member->run_id maps;
//      then the population runs after (nstack..).   <- the chance run-id layout
//   3. RUN: mark every run complete via update_run, in a random "agent completion" order (write burst)
//   4. PROCESS: for EACH member -> a full get_nruns + get_run_status scan over ALL runs (this is what
//      get_failed_run_ids does, once per member) + get_run/get_info over that member's stack run_ids.
// Build + run on Linux (libstdc++); on mac (libc++) it passes.
//   Usage: rs_rwtest [niter] [seed0]
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "RunStorage.h"
using namespace std;

static vector<char> serial_buf(int npar, int nobs) { return vector<char>((npar + nobs) * sizeof(double), 0); }

// one generation's worth of nested-chance queue/run/process; throws on a RunStorage read failure
static void run_once(unsigned seed, int nmember, int nstack_per, int npar, int nobs,
                     const vector<string> &par_names, const vector<string> &obs_names)
{
    std::mt19937 rng(seed);
    RunStorage rs("rs_rwtest.rns");
    rs.reset(par_names, obs_names, "rs_rwtest.rns");          // (1) reinitialize
    vector<double> mp(npar, 1.0);

    // (2) QUEUE: stack runs first, per member, then the population runs
    map<int, vector<int>> member_stack_runs;                 // member -> stack run_ids (mirrors the run-id map)
    for (int m = 0; m < nmember; m++)
    {
        vector<int> rids;
        for (int s = 0; s < nstack_per; s++)
            rids.push_back(rs.add_run(mp, "GEN=0_MEMBER=" + to_string(m), 0.0));
        member_stack_runs[m] = rids;
    }
    vector<int> pop_runs;
    for (int m = 0; m < nmember; m++)
        pop_runs.push_back(rs.add_run(mp, "POP_MEMBER=" + to_string(m), 0.0));

    const int total = nmember * nstack_per + nmember;

    // (3) RUN: store every result via update_run in a random order (stand-in for agent completion order)
    vector<int> all_ids;
    for (int i = 0; i < total; i++) all_ids.push_back(i);
    std::shuffle(all_ids.begin(), all_ids.end(), rng);
    vector<char> sd = serial_buf(npar, nobs);
    for (int rid : all_ids) rs.update_run(rid, sd);

    // (4) PROCESS: per member -> full get_run_status scan (get_failed_run_ids) + per-member reads
    for (int m = 0; m < nmember; m++)
    {
        int n = rs.get_nruns();                              // get_failed_run_ids reads n_runs first
        for (int id = 0; id < n; id++) rs.get_run_status(id);// ...then scans every run's status
        for (int rid : member_stack_runs[m])                 // ...then reads this member's stack runs
        {
            vector<double> pv, ov; string itx; double iv;
            rs.get_run(rid, pv, ov, itx, iv);
            int st; rs.get_info(rid, st, itx, iv);
        }
    }
    rs.free_memory();
}

int main(int argc, char **argv)
{
    const int nmember = 10, nstack_per = 10, npar = 6, nobs = 4;
    const long niter = (argc > 1) ? atol(argv[1]) : 100000L;
    const unsigned seed0 = (argc > 2) ? (unsigned)strtoul(argv[2], nullptr, 10) : 0u;

    vector<string> par_names, obs_names;
    for (int i = 0; i < npar; i++) par_names.push_back("par" + to_string(i));
    for (int i = 0; i < nobs; i++) obs_names.push_back("obs" + to_string(i));

    for (long it = 0; it < niter; it++)
    {
        if (it % 5000 == 0)
            cerr << "\r  " << it << " / " << niter << " gens   " << flush;
        try { run_once((unsigned)(seed0 + it), nmember, nstack_per, npar, nobs, par_names, obs_names); }
        catch (const exception &e)
        {
            cout << "\nREPRODUCED at gen/seed " << (seed0 + it) << ":  " << e.what()
                 << "   (replay: rs_rwtest 1 " << (seed0 + it) << ")" << endl;
            return 1;
        }
    }
    cout << "\nno failure in " << niter << " gens (this build's stdlib tolerates the pattern)" << endl;
    return 0;
}
