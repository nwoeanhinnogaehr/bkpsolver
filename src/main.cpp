#include "comb_solver.hpp"
#include "dcs_solver.hpp"
#include <clipp.h>
#include <iostream>
using namespace std;
using namespace clipp;

int main(int argc, char **argv) {
    int alpha = 2, beta = 2, delta = 30, mu = 250, gamma = 5, omega = 3;
    int lookback = 1;
    bool quiet, use_file, best_prefix = false, weak_lb = false, lb_only = false;
    int num_threads = 4;
    string filename;
    enum class mode {dcs, comb};
    mode selected = mode::comb;
    auto cli = (
        command("dcs").set(selected, mode::dcs)
            % "use algorithm by Federico Della Croce and Rosario Scatamacchia"
            & ( option("-a") & integer("alpha=2", alpha),
                option("-b") & integer("beta=2", beta),
                option("-g") & integer("gamma=5", gamma),
                option("-d") & integer("delta=30", delta),
                option("-m") & integer("mu=250", mu),
                option("-o") & integer("omega=3", omega) ) |
        command("comb").set(selected, mode::comb)
            % "use combinatorial algorithm"
            & ( option("-l") & integer("lookback=1", lookback),
                option("-p").set(best_prefix) % "search for best prefix (try this if number of leaves is large)",
                option("-w").set(weak_lb),
                option("-lb-only").set(lb_only) ),
        option("-q").set(quiet) % "quiet mode: do not log to stderr",
        option("-j") & integer("num_threads=4", num_threads),
        opt_value("file", filename).set(use_file) % "filename (will read from stdin if absent)"
    );
    if(!parse(argc, argv, cli)) {
        cout << make_man_page(cli, argv[0]);
        return 0;
    }

    BKPInstance inst = use_file
        ? BKPInstance::from_file(filename)
        : BKPInstance::from_stdin();

    Comb_BKPSolver<uint16_t> comb_solver_16(inst);
    comb_solver_16.log = !quiet;
    comb_solver_16.lookback = lookback;
    comb_solver_16.best_prefix = best_prefix;
    comb_solver_16.num_threads = num_threads;
    comb_solver_16.weak_lb = weak_lb;
    comb_solver_16.lb_only = lb_only;

    Comb_BKPSolver<uint32_t> comb_solver_32(inst);
    comb_solver_32.log = !quiet;
    comb_solver_32.lookback = lookback;
    comb_solver_32.best_prefix = best_prefix;
    comb_solver_32.num_threads = num_threads;
    comb_solver_32.weak_lb = weak_lb;
    comb_solver_32.lb_only = lb_only;

    DCS_BKPSolver dcs_solver(inst);
    dcs_solver.log = !quiet;
    dcs_solver.alpha = alpha;
    dcs_solver.beta = beta;
    dcs_solver.delta = delta;
    dcs_solver.mu = mu;
    dcs_solver.gamma = gamma;
    dcs_solver.omega = omega;
    dcs_solver.num_threads = num_threads;

    BKPSolution<max_pft_t> sol(inst.n);
    if (selected == mode::comb) {
        try {
            sol = BKPSolution<max_pft_t>(comb_solver_16.solve());
        } catch (ProfitOverflowError) {
            if (!quiet)
                cerr << "profit overflowed 16 bits, retrying with 32 bits" << endl;
            sol = BKPSolution<max_pft_t>(comb_solver_32.solve());
        }
    } else
        sol = dcs_solver.solve();

    if (!quiet && !lb_only) {
        cerr << "upper: ";
        for (size_t i = 0; i < inst.n; i++)
            cerr << (int)sol.up_sol[i];
        cerr << endl << "lower: ";
        for (size_t i = 0; i < inst.n; i++)
            cerr << (int)sol.lo_sol[i];
        cerr << endl;
        cerr << "profit: ";
    }
    cout << sol.pft << endl;
}
