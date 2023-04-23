#include "comb_solver.hpp"
#include "dcs_solver.hpp"
#include <clipp.h>
#include <iostream>
#include <sys/resource.h>
using namespace std;
using namespace clipp;

int main(int argc, char **argv) {
    int alpha = 2, beta = 2, delta = 30, mu = 250, gamma = 5, omega = 3;
    int lookback = 0;
    bool quiet, use_file, best_prefix = false, lb_only = false;
    double density_threshold = 0.02, prefix_density_threshold = 0.02;
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
            & ( (option("-l") & integer("lookback=0", lookback)) % "try deleting items up to this far back from the end of the greedy lower solution",
                option("-p").set(best_prefix) % "search for best prefix (try this if number of leaves is large)",
                (option("-dt") & number("density_threshold=0.02", density_threshold)) % "density threshold for switching to dense DP table",
                (option("-pdt") & number("prefix_density_threshold=0.02", prefix_density_threshold)) % "density threshold for switching to dense prefix DP table",
                option("-lb-only").set(lb_only) % "compute lower bound only"),
        option("-q").set(quiet) % "quiet mode: do not log to stderr",
        (option("-j") & integer("num_threads=4", num_threads)) % "number of threads to use for multithreaded code",
        opt_value("file", filename).set(use_file) % "filename (will read from stdin if absent)"
    );
    if(!parse(argc, argv, cli)) {
        cout << make_man_page(cli, argv[0]);
        return 0;
    }

    cout << "algorithm " << (selected == mode::comb ? "comb" : "dcs") << endl;
    cout << "file " << (use_file ? filename : "-") << endl;
    cout << "threads " << num_threads << endl;
    if (selected == mode::comb) {
        cout << "lookback " << lookback << endl;
        cout << "best_prefix " << best_prefix << endl;
        cout << "lb_only " << lb_only << endl;
        cout << "density_threshold " << density_threshold << endl;
        cout << "prefix_density_threshold " << prefix_density_threshold << endl;
    } else {
        cout << "alpha " << alpha << endl;
        cout << "beta " << beta << endl;
        cout << "gamma " << gamma << endl;
        cout << "delta " << delta << endl;
        cout << "mu " << mu << endl;
        cout << "omega " << omega << endl;
    }

    BKPInstance inst = use_file
        ? BKPInstance::from_file(filename)
        : BKPInstance::from_stdin();

    Comb_BKPSolver<uint16_t> comb_solver_16(inst);
    comb_solver_16.log = !quiet;
    comb_solver_16.lookback = lookback;
    comb_solver_16.best_prefix = best_prefix;
    comb_solver_16.num_threads = num_threads;
    comb_solver_16.lb_only = lb_only;
    comb_solver_16.density_threshold = density_threshold;
    comb_solver_16.prefix_density_threshold = prefix_density_threshold;

    Comb_BKPSolver<uint32_t> comb_solver_32(inst);
    comb_solver_32.log = !quiet;
    comb_solver_32.lookback = lookback;
    comb_solver_32.best_prefix = best_prefix;
    comb_solver_32.num_threads = num_threads;
    comb_solver_32.lb_only = lb_only;
    comb_solver_32.density_threshold = density_threshold;
    comb_solver_32.prefix_density_threshold = prefix_density_threshold;

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
            cout << "pft_bits 16" << endl;
        } catch (ProfitOverflowError) {
            if (!quiet)
                cerr << "profit overflowed 16 bits, retrying with 32 bits" << endl;
            sol = BKPSolution<max_pft_t>(comb_solver_32.solve());
            cout << "pft_bits 32" << endl;
        }
    } else
        sol = dcs_solver.solve();

    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    cout << "memory " << rusage.ru_maxrss/1024.0 << endl;

    if (!lb_only) {
        cout << "upper ";
        int upper_util = 0;
        for (size_t i = 0; i < inst.n; i++) {
            cout << (int)sol.up_sol[i];
            upper_util += sol.up_sol[i] * inst.up_wt[i];
        }
        cout << endl << "lower ";
        int lower_util = 0;
        for (size_t i = 0; i < inst.n; i++) {
            cout << (int)sol.lo_sol[i];
            lower_util += sol.lo_sol[i] * inst.lo_wt[i];
        }
        cout << endl;
        cout << "up_util " << upper_util << endl;
        cout << "lo_util " << lower_util << endl;
        cout << "profit ";
    }
    cout << sol.pft << endl;
}
