#pragma once

#include "bkp_instance.hpp"

#include "knapsacksolver/solution.hpp"
#include "knapsacksolver/algorithms/algorithms.hpp"
#include <gurobi_c++.h>

#include <vector>

using max_pft_t = uint32_t;

template<typename P>
struct BKPSolution {
    P pft;
    std::vector<uint8_t> up_sol, lo_sol;

    BKPSolution(size_t n) : pft(0), up_sol(n), lo_sol(n) {}
    template<typename Q>
    BKPSolution(BKPSolution<Q> q) : pft(q.pft), up_sol(q.up_sol), lo_sol(q.lo_sol) {
        ASSERT(max_pft_t(q.pft) == max_pft_t(pft));
    }
};

template<typename P>
class BKPSolver {
public:
    BKPSolver(BKPInstance inst) : inst(inst), log(true) { }
    virtual ~BKPSolver() {}
    virtual BKPSolution<P> solve() = 0;

    bool log = true;

protected:
    void sort_inst();
    void unsort_sol(BKPSolution<P> &sol);
    int find_crit_lb();
    int find_crit_ub();
    int solve_knapsack(size_t n, std::vector<int> pft, std::vector<int> weight, int cap);

    BKPInstance inst;
    std::vector<size_t> pi;
};

template<typename P>
int BKPSolver<P>::solve_knapsack(size_t n, std::vector<int> pft, std::vector<int> weight, int cap) {
    knapsacksolver::Instance kinst;
    kinst.set_capacity(cap);
    for (size_t i = 0; i < n; i++)
        kinst.add_item(weight[i], pft[i]);
    knapsacksolver::DynamicProgrammingPrimalDualOptionalParameters param;
    param.set_combo();
    return knapsacksolver::dynamic_programming_primal_dual(kinst, param).solution.profit();
}
template<typename P>
void BKPSolver<P>::sort_inst() {
    BKPInstance copy = inst;
    pi.resize(inst.n);
    iota(pi.begin(), pi.end(), 0);
    sort(pi.begin(), pi.end(), [&](int i, int j) {
            // sort first by efficiency, breaking ties by profit
            if (inst.pft[i]*inst.lo_wt[j] == inst.pft[j]*inst.lo_wt[i])
                return inst.pft[i] > inst.pft[j];
            return inst.pft[i]*inst.lo_wt[j] > inst.pft[j]*inst.lo_wt[i];
        });
    for (int i = 0; i < inst.n; i++) {
        inst.up_wt[i] = copy.up_wt[pi[i]];
        inst.lo_wt[i] = copy.lo_wt[pi[i]];
        inst.pft[i] = copy.pft[pi[i]];
    }
}
template<typename P>
void BKPSolver<P>::unsort_sol(BKPSolution<P> &sol) {
    BKPSolution<P> copy = sol;
    for (int i = 0; i < inst.n; i++) {
        sol.lo_sol[pi[i]] = copy.lo_sol[i];
        sol.up_sol[pi[i]] = copy.up_sol[i];
    }
}
template<typename P>
int BKPSolver<P>::find_crit_lb() {
    int w = 0;
    for (int i = 0; i < inst.n; i++) {
        w += inst.lo_wt[i];
        if (w >= inst.lo_cap)
            return i;
    }
    return inst.n;
}
template<typename P>
int BKPSolver<P>::find_crit_ub() {
    int z = solve_knapsack(inst.n, inst.lo_wt, inst.up_wt, inst.up_cap);
    int w = 0;
    for (int i = 0; i < inst.n; i++) {
        w += inst.lo_wt[i];
        if (w >= inst.lo_cap + z)
            return i;
    }
    return inst.n;
}
