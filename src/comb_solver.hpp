#pragma once

#include "bkp_instance.hpp"
#include "bkp_solver.hpp"

#include "knapsacksolver/solution.hpp"
#include "knapsacksolver/algorithms/minknap.hpp"
#include <gurobi_c++.h>
#include <omp.h>

#include <vector>
#include <iostream>
#include <chrono>

struct ProfitOverflowError {};

template<typename P>
class Comb_BKPSolver : public BKPSolver<P> {
public:
    Comb_BKPSolver(BKPInstance inst) : BKPSolver<P>(inst), grb_env(true), ub(inst.n) { }

    BKPSolution<P> solve() override;
    virtual ~Comb_BKPSolver() override {
        if (lb_table) free(lb_table);
        if (lb_pre_tab) free(lb_pre_tab);
    }

    int lookback = 1;
    bool best_prefix = false;
    int num_threads = 4;
    bool weak_lb = false;
    bool lb_only = false;
    using BKPSolver<P>::log;

private:
    void solve_rec(size_t i, int up_cap, int lo_cap, int pft);
    BKPSolution<P> solve_lower();
    bool lower_bound_test(size_t i, int up_cap, int lo_cap, int pft);
    void compute_lower_bound();
    void compute_lower_bound_weak();
    P &lower_bound(size_t i, int up_cap, int lo_cap);
    BKPSolution<max_pft_t> greedy_ub();
    double dcs_lb();
    double dcs_lb(int c);

    GRBEnv grb_env;
    P *lb_table = nullptr;
    P *lb_pre_tab = nullptr;
    BKPSolution<P> ub;
    std::vector<uint8_t> cur_up_sol;
    std::vector<size_t> lo_greedy_items;
    std::vector<int> min_upper, min_lower;
    uint64_t nodes = 0, leaves = 0, max_depth = 0, sum_depth = 0;
    std::chrono::system_clock::time_point opt_time;

    using BKPSolver<P>::inst;
};

template<typename P>
BKPSolution<P> Comb_BKPSolver<P>::solve() {
    using namespace std::chrono;

    auto init_time = high_resolution_clock::now();
    
    grb_env.set("OutputFlag", "0");
    grb_env.start();

    BKPSolver<P>::sort_inst();

    if (!lb_only) {
        if (log) std::cerr << "initial bound test... ";
        BKPSolution<max_pft_t> initial_ub = greedy_ub();
        if (initial_ub.pft >= (1ull<<(8*sizeof(P))))
            throw ProfitOverflowError{};
        std::cout << "greedy_ub " << initial_ub.pft << std::endl;
        P lb = dcs_lb();
        std::cout << "dcs_lb " << lb << std::endl;
        auto initial_bound_end_time = high_resolution_clock::now();
        std::cout << "initial_bound_test_time " << duration<double, std::milli>(initial_bound_end_time - init_time).count() << std::endl;
        if (initial_ub.pft <= lb) {
            if (log) std::cerr << "done, returning early" << std::endl;
            std::cout << "total_time " << duration<double, std::milli>(initial_bound_end_time - init_time).count() << std::endl;
            return initial_ub;
        }
        ub = initial_ub;
    }

    auto lb_start_time = high_resolution_clock::now();
    if (log) std::cerr << "computing lower bounds... ";
    if (weak_lb) compute_lower_bound_weak();
    else compute_lower_bound();
    if (lb_only) {
        ub.pft = lower_bound(0, inst.up_cap, inst.lo_cap);
        return ub;
    }
    std::cout << "dp_lb " << lower_bound(0, inst.up_cap, inst.lo_cap) << std::endl;
    auto lb_end_time = high_resolution_clock::now();
    std::cout << "dp_lb_time " << duration<double, std::milli>(lb_end_time - lb_start_time).count() << std::endl;

    opt_time = lb_end_time;
    if (log) std::cerr << "starting search... " << std::endl;
    cur_up_sol.resize(inst.n+1);
    fill(cur_up_sol.begin(), cur_up_sol.end(), 1);
    solve_rec(0, inst.up_cap, inst.lo_cap, 0);
    std::cout << "nodes " << nodes << std::endl;
    std::cout << "leaves " << leaves << std::endl;

    BKPSolver<P>::unsort_sol(ub);

    auto bnb_end_time = high_resolution_clock::now();
    std::cout << "bnb_time " << duration<double, std::milli>(bnb_end_time - lb_end_time).count() << std::endl;
    std::cout << "opt_time " << duration<double, std::milli>(opt_time - lb_end_time).count() << std::endl;
    std::cout << "proof_time " << duration<double, std::milli>(bnb_end_time - opt_time).count() << std::endl;
    std::cout << "total_time " << duration<double, std::milli>(bnb_end_time - init_time).count() << std::endl;

    std::cout << "max_depth " << max_depth << std::endl;
    std::cout << "avg_depth " << sum_depth/(double)nodes << std::endl;
 
    return ub;
}

template<typename P>
void Comb_BKPSolver<P>::solve_rec(size_t i, int up_cap, int lo_cap, int pft) {
    max_depth = std::max(max_depth, i);
    sum_depth += i;

    if (i == inst.n) {
        leaves++;
        BKPSolution<P> resp = solve_lower();
        if (resp.pft < ub.pft) {
            ub = resp;
            if (log) std::cerr << "new upper bound: " << ub.pft << std::endl;
            opt_time = std::chrono::high_resolution_clock::now();
        }
        return;
    }

    nodes++;

    if (log) if (nodes % (1<<20) == 0) {
        std::cerr << (nodes >> 20) << "M nodes and " << leaves << " leaves searched" << std::endl;
        std::cerr << "current node: ";
        double p = 1;
        for (size_t j = 0; j < std::min(i,size_t(64)); j++) {
            std::cerr << (int)cur_up_sol[j];
            if (cur_up_sol[j])
                p -= pow(2.0, -int(j)-1);
        }
        for (int j = i; j <= inst.n; j++)
            p -= pow(2.0, -int(j)-1);
        if (i > 64)
            std::cerr << "... depth=" << i;
        std::cerr << std::endl;
        std::cerr << "progress: " << p*100 << "%" << std::endl;
    }

    if (lower_bound_test(i, up_cap, lo_cap, pft))
        return;

    if (inst.up_wt[i] <= up_cap)
            solve_rec(i+1, up_cap-inst.up_wt[i], lo_cap, pft);

    cur_up_sol[i] = 0;
    if (inst.lo_wt[i] <= lo_cap) {
        lo_greedy_items.push_back(i);
        solve_rec(i+1, up_cap, lo_cap-inst.lo_wt[i], pft+inst.pft[i]);
        lo_greedy_items.pop_back();
    }
    else solve_rec(i+1, up_cap, lo_cap, pft);
    cur_up_sol[i] = 1;
}

template<typename P>
BKPSolution<P> Comb_BKPSolver<P>::solve_lower() {
    knapsacksolver::Instance finst;
    finst.set_capacity(inst.lo_cap);
    for (int i = 0; i < inst.n; i++)
        if (!cur_up_sol[i])
            finst.add_item(inst.lo_wt[i], inst.pft[i]);
    knapsacksolver::MinknapOptionalParameters param;
    param.set_combo();
    auto foutput = knapsacksolver::minknap(finst, param);
    auto fsol = foutput.solution;
    BKPSolution<P> bkpsol(inst.n);
    bkpsol.pft = fsol.profit();
    bkpsol.up_sol = cur_up_sol;
    for (int i = 0, j = 0; i < inst.n; i++)
        if (!cur_up_sol[i]) {
            if (fsol.contains_idx(j)) bkpsol.lo_sol[i] = 1;
            else bkpsol.lo_sol[i] = 0;
            j++;
        }
    return bkpsol;
}

template<typename P>
bool Comb_BKPSolver<P>::lower_bound_test(size_t i, int up_cap, int lo_cap, int pft) {
    if (pft >= ub.pft)
        return true;
    if (lower_bound(i, up_cap, lo_cap)+pft >= ub.pft)
        return true;
    int limit = std::max(0, int(lo_greedy_items.size())-lookback);
    for (int j = int(lo_greedy_items.size())-1; j >= limit; j--) {
        size_t k = lo_greedy_items[j];
        if (lower_bound(i, up_cap, lo_cap+inst.lo_wt[k])+pft-inst.pft[k] >= ub.pft)
            return true;
    }
    if (best_prefix && i > 0) {
        int next = -1;
        for (int j = i-2; j >= 0; j--)
            if (cur_up_sol[j] == 0) {
                next = j;
                break;
            }
        if (cur_up_sol[i-1] == 0) {
            if (next < 0) {
                for (int k = 0; k < inst.lo_wt[i-1]; k++)
                    lb_pre_tab[(i-1)*(inst.lo_cap+1)+k] = 0;
                for (int k = inst.lo_wt[i-1]; k <= inst.lo_cap; k++)
                    lb_pre_tab[(i-1)*(inst.lo_cap+1)+k] = inst.pft[i-1];
            } else {
                for (int k = 0; k < inst.lo_wt[i-1]; k++)
                    lb_pre_tab[(i-1)*(inst.lo_cap+1)+k] = lb_pre_tab[next*(inst.lo_cap+1)+k];
                for (int k = inst.lo_wt[i-1]; k <= inst.lo_cap; k++)
                    lb_pre_tab[(i-1)*(inst.lo_cap+1)+k] = std::max(
                            lb_pre_tab[next*(inst.lo_cap+1)+k],
                            P(lb_pre_tab[next*(inst.lo_cap+1)+k-inst.lo_wt[i-1]]+inst.pft[i-1]));
            }
            next = i-1;
        }
        if (next >= 0) {
            /* Normally I wouldn't go to such extremes to optimize something,
               but this is the hottest loop in the program for many hard instances.
               Original code:

            for (int k = min_lower[i]; k <= inst.lo_cap; k++)
                if (lb_pre_tab[next*(inst.lo_cap+1) + inst.lo_cap - k]
                        + lower_bound(i, up_cap, k) >= ub.pft)
                    return true;
            return false;

            The below SIMD-friendly implementation gives a speedup of up to 8x. */
            if (lb_pre_tab[next*(inst.lo_cap+1)+inst.lo_cap] >= ub.pft)
                return true;
            const int block_size = 256;
            int k = 0;
            P best = 0;
            for (; k <= inst.lo_cap-min_lower[i]-block_size+1; k += block_size) {
                for (int j = k; j < k+block_size; j++) {
                    best = std::max(best, P(lb_pre_tab[next*(inst.lo_cap+1)+j]
                            + lower_bound(i, up_cap, inst.lo_cap-j)));
                }
                if (best >= ub.pft)
                    return true;
            }
            for (int j = k; j <= inst.lo_cap-min_lower[i]; j++) {
                best = std::max(best, P(lb_pre_tab[next*(inst.lo_cap+1)+j]
                        + lower_bound(i, up_cap, inst.lo_cap-j)));
            }
            if (best >= ub.pft)
                return true;
        }
    }

    return false;
}

template<typename P>
void Comb_BKPSolver<P>::compute_lower_bound() {
    lb_pre_tab = (P*) malloc((inst.n+1)*(inst.lo_cap+1)*sizeof(P));
    ASSERT(lb_pre_tab);
    lb_table = (P*) malloc((inst.n+1)*(inst.up_cap+1)*(inst.lo_cap+1)*sizeof(P));
    ASSERT(lb_table);
    min_upper.resize(inst.n);
    min_lower.resize(inst.n);
    min_upper[0] = inst.up_cap;
    min_lower[0] = inst.lo_cap;
    for (size_t i = 1; i < inst.n; i++) {
        min_upper[i] = std::max(0, min_upper[i-1] - inst.up_wt[i-1]);
        min_lower[i] = std::max(0, min_lower[i-1] - inst.lo_wt[i-1]);
    }
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int ru = 0; ru <= inst.up_cap; ru++)
            for (int rl = 0; rl <= inst.lo_cap; rl++)
                lower_bound(inst.n, ru, rl) = 0;
        for (int i = inst.n-1; i >= 0; i--) {
            #pragma omp barrier
            #pragma omp for nowait
            for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
                for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                    lower_bound(i, ru, rl) = std::min(lower_bound(i+1, ru-inst.up_wt[i], rl),
                            std::max(P(lower_bound(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i]),
                                lower_bound(i+1, ru, rl)));

            if (min_lower[i] <= inst.lo_wt[i] && min_upper[i] <= inst.up_wt[i])
                #pragma omp for nowait
                for (int ru = 0; ru < inst.up_wt[i]; ru++)
                    for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                        lower_bound(i, ru, rl) = lower_bound(i+1, ru, rl);

            if (min_lower[i] <= inst.lo_wt[i])
                #pragma omp for nowait
                for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
                    for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                        lower_bound(i, ru, rl) = std::min(lower_bound(i+1, ru-inst.up_wt[i], rl),
                                lower_bound(i+1, ru, rl));

            if (min_upper[i] <= inst.up_wt[i])
                #pragma omp for nowait
                for (int ru = 0; ru < inst.up_wt[i]; ru++)
                    for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                        lower_bound(i, ru, rl) = std::max(
                                P(lower_bound(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i]),
                                lower_bound(i+1, ru, rl));
        }
    }
}

template<typename P>
void Comb_BKPSolver<P>::compute_lower_bound_weak() {
    lb_pre_tab = (P*) malloc((inst.n+1)*(inst.lo_cap+1)*sizeof(P));
    ASSERT(lb_pre_tab);
    lb_table = (P*) malloc((inst.n+1)*(inst.up_cap+1)*(inst.lo_cap+1)*sizeof(P));
    ASSERT(lb_table);
    min_upper.resize(inst.n);
    min_lower.resize(inst.n);
    min_upper[0] = inst.up_cap;
    min_lower[0] = inst.lo_cap;
    for (size_t i = 1; i < inst.n; i++) {
        min_upper[i] = std::max(0, min_upper[i-1] - inst.up_wt[i-1]);
        min_lower[i] = std::max(0, min_lower[i-1] - inst.lo_wt[i-1]);
    }
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int ru = 0; ru <= inst.up_cap; ru++)
            for (int rl = 0; rl <= inst.lo_cap; rl++)
                lower_bound(inst.n, ru, rl) = 0;
        for (int i = inst.n-1; i >= 0; i--) {
            #pragma omp barrier
            #pragma omp for nowait
            for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
                for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                    lower_bound(i, ru, rl) = std::min(lower_bound(i+1, ru-inst.up_wt[i], rl),
                            P(lower_bound(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i]));

            if (min_lower[i] <= inst.lo_wt[i] && min_upper[i] <= inst.up_wt[i])
                #pragma omp for nowait
                for (int ru = 0; ru < inst.up_wt[i]; ru++)
                    for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                        lower_bound(i, ru, rl) = lower_bound(i+1, ru, rl);

            if (min_lower[i] <= inst.lo_wt[i])
                #pragma omp for nowait
                for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
                    for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                        lower_bound(i, ru, rl) = std::min(lower_bound(i+1, ru-inst.up_wt[i], rl),
                                lower_bound(i+1, ru, rl));

            if (min_upper[i] <= inst.up_wt[i])
                #pragma omp for nowait
                for (int ru = 0; ru < inst.up_wt[i]; ru++)
                    for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                        lower_bound(i, ru, rl) = lower_bound(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i];
        }
    }
}

template<typename P>
P &Comb_BKPSolver<P>::lower_bound(size_t i, int up_cap, int lo_cap) {
    return lb_table[i*(inst.up_cap+1)*(inst.lo_cap+1) + up_cap*(inst.lo_cap+1) + lo_cap];
}

template<typename P>
BKPSolution<max_pft_t> Comb_BKPSolver<P>::greedy_ub() {
    BKPSolution<max_pft_t> bkpsol(inst.n);

    knapsacksolver::Instance kpinst;
    kpinst.set_capacity(inst.up_cap);
    for (int i = 0; i < inst.n; i++)
        kpinst.add_item(inst.up_wt[i], inst.pft[i]);
    knapsacksolver::MinknapOptionalParameters param;
    param.set_combo();
    auto output = knapsacksolver::minknap(kpinst, param);
    auto sol = output.solution;

    knapsacksolver::Instance finst;
    finst.set_capacity(inst.lo_cap);
    for (int i = 0; i < inst.n; i++) {
        if (!sol.contains_idx(i)) {
            finst.add_item(inst.lo_wt[i], inst.pft[i]);
            bkpsol.up_sol[i] = 0;
        } else {
            finst.add_item(inst.lo_wt[i], 0);
            bkpsol.up_sol[i] = 1;
        }
    }
    auto foutput = knapsacksolver::minknap(finst, param);
    auto fsol = foutput.solution;
    for (int i = 0; i < inst.n; i++) {
        if (sol.contains_idx(i)) continue;
        if (fsol.contains_idx(i)) bkpsol.lo_sol[i] = 1;
        else bkpsol.lo_sol[i] = 0;
    }
    bkpsol.pft = fsol.profit();
    return bkpsol;
}

template<typename P>
double Comb_BKPSolver<P>::dcs_lb() {
    grb_env.set("Threads", "1");

    int clb = BKPSolver<P>::find_crit_lb();
    int cub = BKPSolver<P>::find_crit_ub();
    double best = INFINITY;
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int c = clb; c <= std::min(int(inst.n)-1, cub); c++) {
        double lb = dcs_lb(c);
        #pragma omp critical
        best = std::min(best, lb);
    }
    return best;
}

template<typename P>
double Comb_BKPSolver<P>::dcs_lb(int c) {
    GRBModel model(grb_env);

    std::vector<GRBVar> vars;
    GRBLinExpr obj = 0.0;
    for (int i = 0; i < c; i++) {
        GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x" + std::to_string(i));
        vars.push_back(var);
        obj += inst.pft[i] - var*inst.pft[i];
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBLinExpr up_kp = 0.0;
    for (int i = 0; i < c; i++)
        up_kp += vars[i]*inst.up_wt[i];
    model.addConstr(up_kp <= inst.up_cap, "up_kp");
    GRBLinExpr lo_kp = 0.0;
    for (int i = 0; i < c; i++)
        lo_kp += inst.lo_wt[i] - vars[i]*inst.lo_wt[i];
    model.addConstr(lo_kp <= inst.lo_cap, "lo_kp");
    model.addConstr(lo_kp >= inst.lo_cap - inst.lo_wt[c] + 1, "lo_crit");

    model.optimize();
    if (model.get(GRB_IntAttr_SolCount) == 0)
        return INFINITY;
    return model.get(GRB_DoubleAttr_ObjVal);
}
