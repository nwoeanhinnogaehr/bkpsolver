#pragma once

#include "bkp_instance.hpp"
#include "bkp_solver.hpp"

#include "knapsacksolver/solution.hpp"
#include "knapsacksolver/algorithms/algorithms.hpp"
#include <gurobi_c++.h>
#include <omp.h>

#include <algorithm>
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
        if (d_tab) free(d_tab);
        if (d_pre_tab) free(d_pre_tab);
    }

    int lookback = 0;
    bool best_prefix = true;
    int num_threads = 4;
    bool lb_only = false;
    double density_threshold = 0.02, prefix_density_threshold = 0.02;
    using BKPSolver<P>::log;

private:
    // sparse tables
    struct __attribute__((packed)) LowerTableEntry {
        int wt;
        P pft;
        bool operator==(const LowerTableEntry& other) {
            return wt == other.wt && pft == other.pft;
        }
    };
    struct __attribute__((packed)) UpperTableEntry {
        int wt;
        std::vector<LowerTableEntry> L;
    };

    void solve_rec(size_t i, int up_cap, int lo_cap, int pft);
    BKPSolution<P> solve_lower();
    bool lower_bound_test_dense(size_t i, int up_cap, int lo_cap, int pft);
    bool lower_bound_test(size_t i, int up_cap, int lo_cap, int pft);
    void init_capacity_bounds();
    void compute_lower_bound_sparse();
    void compute_lower_bound_dense(size_t start_at);
    P &d_tab_val(size_t i, int up_cap, int lo_cap); // dense DP table value O(1)
    P s_tab_val(size_t i, int up_cap, int lo_cap); // sparse DP table value O(log(up_cap) log(lo_cap))
    P lower_bound(size_t i, int up_cap, int lo_cap);
    std::vector<UpperTableEntry> &get_item_table_sparse(size_t i);
    UpperTableEntry &get_upper_entry_sparse(size_t i, int up_cap);
    BKPSolution<max_pft_t> greedy_ub();
    double dcs_lb();
    double dcs_lb(int c);
    void dump_lower_bound();

    GRBEnv grb_env;
    BKPSolution<P> ub;
    std::vector<uint8_t> cur_up_sol;
    std::vector<size_t> lo_greedy_items;
    int n_not_interdicted = 0;
    std::vector<int> min_upper, min_lower;
    uint64_t nodes = 0, leaves = 0, max_depth = 0, sum_depth = 0;
    std::chrono::system_clock::time_point opt_time;

    using BKPSolver<P>::inst;

    // sparse tables
    std::vector<std::vector<UpperTableEntry>> s_tab;
    std::vector<std::vector<LowerTableEntry>> s_pre_tab;

     // dense tables
    P *d_tab = nullptr;
    size_t dense_starts_at;
    P *d_pre_tab = nullptr;
    size_t d_pre_tab_idx;

    // given an item, indicates that it uses dense (false) or sparse (true) DP
    std::vector<bool> is_sparse;
    std::vector<bool> pre_is_sparse;
};

template<typename P>
BKPSolution<P> Comb_BKPSolver<P>::solve() {
    using namespace std::chrono;

    auto init_time = high_resolution_clock::now();
    
    grb_env.set("OutputFlag", "0");
    grb_env.start();

    BKPSolver<P>::sort_inst();

    if (!lb_only) {
        if (log) std::cerr << "initial bound test... " << std::endl;
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
    is_sparse.resize(inst.n);
    s_tab.reserve(inst.n+1);
    init_capacity_bounds();
    compute_lower_bound_sparse();
    // dump_lower_bound();
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
    pre_is_sparse.push_back(1); // start out sparse
    d_pre_tab_idx = 0;
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

    if (log) if (nodes % 1000000 == 0) {
        std::cerr << nodes/1000000 << "M nodes and " << leaves << " leaves searched" << std::endl;
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

    // try taking the item at the upper level
    if (inst.up_wt[i] <= up_cap)
        solve_rec(i+1, up_cap-inst.up_wt[i], lo_cap, pft);

    cur_up_sol[i] = 0;
    n_not_interdicted++;
    if (inst.lo_wt[i] <= lo_cap) {
        lo_greedy_items.push_back(i);
        // try taking the item at the lower level
        solve_rec(i+1, up_cap, lo_cap-inst.lo_wt[i], pft+inst.pft[i]);
        lo_greedy_items.pop_back();
    }
    // try ignoring this item
    else solve_rec(i+1, up_cap, lo_cap, pft);
    cur_up_sol[i] = 1;
    n_not_interdicted--;
}

template<typename P>
BKPSolution<P> Comb_BKPSolver<P>::solve_lower() {
    knapsacksolver::Instance finst;
    finst.set_capacity(inst.lo_cap);
    for (int i = 0; i < inst.n; i++)
        if (!cur_up_sol[i])
            finst.add_item(inst.lo_wt[i], inst.pft[i]);
    knapsacksolver::DynamicProgrammingPrimalDualOptionalParameters param;
    param.set_combo();
    auto foutput = knapsacksolver::dynamic_programming_primal_dual(finst, param);
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
        // delete table rows which are from old branches of the search tree
        while (s_pre_tab.size() + d_pre_tab_idx > n_not_interdicted - !cur_up_sol[i-1]) {
            ASSERT(pre_is_sparse.size() >= 2);
            if (pre_is_sparse[pre_is_sparse.size()-2]) {
                ASSERT(!s_pre_tab.empty());
                s_pre_tab.pop_back();
            } else {
                ASSERT(d_pre_tab_idx);
                d_pre_tab_idx--;
            }
            pre_is_sparse.pop_back();

            // if we backtrack once after switching to dense, we could be left with an
            // uninitialized dense table. make sure this doesn't happen.
            if (d_pre_tab_idx == 0)
                pre_is_sparse.back() = 1;
        }

        if (cur_up_sol[i-1] == 0) { // if the previous item wasn't interdicted
            if (pre_is_sparse.back()) {
                if (s_pre_tab.empty()) {
                    // initialize the table with this first item
                    s_pre_tab.push_back({LowerTableEntry{0, 0}, LowerTableEntry{inst.lo_wt[i-1], P(inst.pft[i-1])}});
                } else {
                    ASSERT(pre_is_sparse.size() >= 2 && pre_is_sparse[pre_is_sparse.size()-2]);

                    s_pre_tab.push_back(s_pre_tab.back()); // clone previous row

                    auto &knap = s_pre_tab.back();
                    std::vector<LowerTableEntry> new_knap;
                    new_knap.reserve(knap.size() * 2);
                    for (int j = knap.size() - 1; j >= 0; --j) {
                        LowerTableEntry it = knap[j];
                        int wt = it.wt + inst.lo_wt[i-1];
                        if (wt <= inst.lo_cap) { // if adding this item fits within capacity
                            P pft = it.pft + inst.pft[i-1];
                            int k = knap.size()-1;

                            // copy the undominated items to output (these can never be dominated in the future)
                            while (knap[k].pft >= pft) k--;
                            std::reverse_copy(knap.begin() + k + 1, knap.end(), std::back_inserter(new_knap));

                            // erase the dominated items
                            while (knap[k].wt >= wt) k--;
                            knap.resize(k + 1);

                            // add the new item (as long as it's an improvement)
                            if (new_knap.empty() || new_knap.back().wt > wt)
                                knap.push_back(LowerTableEntry{wt, pft});
                        }
                        // make sure we didn't delete the first item
                        ASSERT(knap[0].wt == 0 && knap[0].pft == 0);
                    }
                    // copy any remaining items
                    std::reverse_copy(knap.begin(), knap.end(), std::back_inserter(new_knap));
                    // new_knap is reversed here
                    ASSERT(new_knap.back().wt == 0 && new_knap.back().pft == 0);

                    // TODO: is there any need to do this, or is it already clean?
                    std::vector<LowerTableEntry> cleaned;
                    cleaned.reserve(new_knap.size());
                    auto last = std::prev(new_knap.end());
                    for (auto it = new_knap.begin(); it != last; ++it)
                        if (std::next(it)->pft != it->pft)
                            cleaned.push_back(*it);
                    cleaned.push_back(new_knap.back());
                    knap.clear();
                    for (auto it = std::prev(cleaned.end()); it != cleaned.begin(); --it)
                        if (std::prev(it)->wt != it->wt)
                            knap.push_back(*it);
                    knap.push_back(cleaned.front());
                }

                double density = s_pre_tab.back().size()/double(inst.lo_cap);
                if (density >= prefix_density_threshold) {
                    pre_is_sparse.back() = 0;
                    pre_is_sparse.push_back(0);
                    // convert sparse -> dense
                    if (!d_pre_tab) {
                        d_pre_tab = (P*) malloc((inst.n+1)*uint64_t(inst.lo_cap+1)*sizeof(P));
                        ASSERT(d_pre_tab && "couldn't allocate memory");
                    }
                    P *dknap = &d_pre_tab[++d_pre_tab_idx*(inst.lo_cap+1)];
                    auto &sknap = s_pre_tab.back();
                    for (auto it = sknap.begin(); it != std::prev(sknap.end()); ++it) {
                        for (int j = it->wt; j < std::next(it)->wt; j++)
                            dknap[j] = it->pft;
                    }
                    ASSERT(sknap.size() >= 1);
                    for (int j = std::prev(sknap.end())->wt; j <= inst.lo_cap; j++)
                        dknap[j] = std::prev(sknap.end())->pft;
                    s_pre_tab.pop_back();
                } else {
                    pre_is_sparse.push_back(1); // continue with sparse
                }
            } else { // use dense table
                if (!d_pre_tab) {
                    d_pre_tab = (P*) malloc((inst.n+1)*uint64_t(inst.lo_cap+1)*sizeof(P));
                    ASSERT(d_pre_tab && "couldn't allocate memory");
                }
                P *knap = &d_pre_tab[++d_pre_tab_idx*(inst.lo_cap+1)];
                if (d_pre_tab_idx == 1) { // initialize table
                    ASSERT(s_pre_tab.size() == 0); // we shouldn't have any table at all if we get here
                    for (int k = 0; k < inst.lo_wt[i-1]; k++)
                        knap[k] = 0;
                    for (int k = inst.lo_wt[i-1]; k <= inst.lo_cap; k++)
                        knap[k] = inst.pft[i-1];
                } else { // update table
                    P *prev = &d_pre_tab[(d_pre_tab_idx-1)*(inst.lo_cap+1)];
                    for (int k = 0; k < inst.lo_wt[i-1]; k++)
                        knap[k] = prev[k];
                    for (int k = inst.lo_wt[i-1]; k <= inst.lo_cap; k++)
                        knap[k] = std::max(prev[k], P(prev[k-inst.lo_wt[i-1]]+inst.pft[i-1]));
                }
                pre_is_sparse.push_back(0); // continue with dense
            }
        }
        if (pre_is_sparse.size() >= 2 && pre_is_sparse[pre_is_sparse.size()-2]) {
            if (!s_pre_tab.empty()) {
                auto &knap = s_pre_tab.back();
                if (is_sparse[i]) {
                    // sparse DP table, sparse prefix DP table
                    UpperTableEntry &uit = get_upper_entry_sparse(i, up_cap);
                    auto e2 = std::prev(uit.L.end());
                    for (LowerTableEntry e1 : knap) {
                        while (e1.wt + e2->wt > inst.lo_cap) --e2;
                        if (e1.pft + e2->pft >= ub.pft)
                            return true;
                    }
                } else {
                    // dense DP table, sparse prefix DP table
                    for (LowerTableEntry e : knap)
                        if (d_tab_val(i, up_cap, inst.lo_cap - e.wt) + e.pft >= ub.pft)
                            return true;
                }
            }
        } else { // dense prefix table
            if (d_pre_tab_idx) {
                P *knap = &d_pre_tab[d_pre_tab_idx*(inst.lo_cap+1)];
                if (is_sparse[i]) {
                    // sparse DP table, dense prefix DP table
                    if (knap[inst.lo_cap] >= ub.pft)
                        return true;
                    const int block_size = 256;
                    int k = 0;
                    P best = 0;
                    UpperTableEntry &uit = get_upper_entry_sparse(i, up_cap);
                    for (; k < int(uit.L.size())-block_size; k += block_size) {
                        for (int j = k; j < k+block_size; j++)
                            best = std::max(best, P(knap[inst.lo_cap-uit.L[j].wt] + uit.L[j].pft));
                        if (best >= ub.pft)
                            return true;
                    }
                    for (int j = k; j < uit.L.size(); j++)
                        best = std::max(best, P(knap[inst.lo_cap-uit.L[j].wt] + uit.L[j].pft));
                    if (best >= ub.pft)
                        return true;
                } else {
                    // dense DP table, dense prefix DP table
                    if (knap[inst.lo_cap] >= ub.pft)
                        return true;
                    const int block_size = 256;
                    int k = 0;
                    P best = 0;
                    for (; k < inst.lo_cap-block_size-min_lower[i]; k += block_size) {
                        for (int j = k; j < k+block_size; j++)
                            best = std::max(best, P(knap[j]
                                    + d_tab_val(i, up_cap, inst.lo_cap-j)));
                        if (best >= ub.pft)
                            return true;
                    }
                    for (int j = k; j <= inst.lo_cap-min_lower[i]; j++)
                        best = std::max(best, P(knap[j]
                                + d_tab_val(i, up_cap, inst.lo_cap-j)));
                    if (best >= ub.pft)
                        return true;
                }
            }
        }
    }

    return false;
}

// calculate a lower bound on the upper/lower level capacity that will be
// needed by branch and bound
template<typename P>
void Comb_BKPSolver<P>::init_capacity_bounds() {
    min_upper.resize(inst.n);
    min_lower.resize(inst.n);
    min_upper[0] = inst.up_cap;
    min_lower[0] = inst.lo_cap;
    for (size_t i = 1; i < inst.n; i++) {
        min_upper[i] = std::max(0, min_upper[i-1] - inst.up_wt[i-1]);
        min_lower[i] = std::max(0, min_lower[i-1] - inst.lo_wt[i-1]);
    }
}

template<typename P>
void Comb_BKPSolver<P>::compute_lower_bound_sparse() {
    // initialize the table with the last item
    s_tab.resize(inst.n);
    s_tab.back().push_back(UpperTableEntry{0, {LowerTableEntry{0,0}, LowerTableEntry{inst.lo_wt[inst.n-1], P(inst.pft[inst.n-1])}}});
    s_tab.back().push_back(UpperTableEntry{inst.up_wt[inst.n-1], {LowerTableEntry{0, 0}}});
    is_sparse[inst.n-1] = 1;

    // loop over all remaining items
    for (int i = inst.n-2; i >= 0; i--) {
        is_sparse[i] = 1;
        // add the new item to the table
        std::vector<UpperTableEntry> prev = s_tab[i+1];
        s_tab[i] = prev;
        auto &cur = s_tab[i];

        // update the columns
        for (auto col = cur.rbegin(); col != cur.rend(); col++) {
            // run lower level knapsack update on column
            auto &knap = col->L;
            std::vector<LowerTableEntry> new_knap;
            new_knap.reserve(knap.size() * 2);
            for (int j = knap.size() - 1; j >= 0; --j) {
                LowerTableEntry it = knap[j];
                int wt = it.wt + inst.lo_wt[i];
                if (wt <= inst.lo_cap) { // if adding this item fits within capacity
                    P pft = it.pft + inst.pft[i];
                    size_t k = knap.size()-1;

                    // copy the undominated items to output (these can never be dominated in the future)
                    while (knap[k].pft >= pft) k--;
                    std::reverse_copy(knap.begin() + k + 1, knap.end(), std::back_inserter(new_knap));

                    // erase the dominated items
                    while (knap[k].wt >= wt) k--;
                    knap.resize(k + 1);

                    // add the new item (as long as it's an improvement)
                    if (new_knap.empty() || new_knap.back().wt > wt)
                        knap.push_back(LowerTableEntry{wt, pft});
                }
                // make sure we didn't delete the first item
                ASSERT(knap[0].wt == 0 && knap[0].pft == 0);
            }
            // copy any remaining items
            std::reverse_copy(knap.begin(), knap.end(), std::back_inserter(new_knap));
            // new_knap is reversed here
            ASSERT(new_knap.back().wt == 0 && new_knap.back().pft == 0);

            // copy it into knap
            knap = new_knap;
        }


        // insert the new columns (cloned)
        size_t ncols_before = cur.size();
        cur.reserve(2 * cur.size());
        auto insert_at = cur.begin();
        auto end = cur.end();
        for (auto it = cur.begin(); it != end; it++) {
            if (it->wt + inst.up_wt[i] <= inst.up_cap) { // if we have enough capacity
                // find the position where the new column gets inserted
                while (insert_at != end && insert_at->wt <= it->wt + inst.up_wt[i])
                    ++insert_at;

                // if we already have this column, ignore it
                if (std::prev(insert_at)->wt == it->wt + inst.up_wt[i]) 
                    continue;

                // add the column
                UpperTableEntry new_col = *std::prev(insert_at);
                new_col.wt = it->wt + inst.up_wt[i];
                cur.push_back(new_col);
            }
        }
        // merge old and new columns
        std::inplace_merge(cur.begin(), cur.begin()+ncols_before, cur.end(),
            [](const UpperTableEntry& a, const UpperTableEntry& b) {
                return a.wt < b.wt;
            }
        );

        auto old_col = std::prev(prev.end());
        for (auto col = cur.rbegin(); col != cur.rend(); col++) {
            auto &knap = col->L;
            // run upper level knapsack update on column, i.e.,
            // replace column entires with the min of that or column we get from interdicting
            std::vector<LowerTableEntry> new_knap;
            new_knap.reserve(2 * knap.size());
            std::vector<LowerTableEntry> *temp = &knap;
            if (col->wt - inst.up_wt[i] >= 0) { // if we have enough interdiction budget
                // find the interdiction column
                while (old_col != prev.begin() && old_col->wt > col->wt - inst.up_wt[i])
                    old_col--;

                auto ptr = knap.begin();
                auto intd_ptr = std::prev(old_col->L.end());
                while (true) {
                    // update profit
                    ptr->pft = std::min(intd_ptr->pft, ptr->pft);
                    // if the interdiction column has more entires, add them
                    if (intd_ptr->wt > ptr->wt) {
                        new_knap.push_back(LowerTableEntry{intd_ptr->wt, ptr->pft});
                        intd_ptr--;
                    } else {
                        // otherwise, just copy the entry to the output
                        new_knap.push_back(*ptr);
                        ptr++;
                        if (ptr == knap.end())
                            break;
                    }
                }
                temp = &new_knap;
            }
            // merge adjacent equal regions
            std::vector<LowerTableEntry> cleaned;
            cleaned.reserve(temp->size());
            auto last = std::prev(temp->end());
            for (auto it = temp->begin(); it != last; ++it)
                if (std::next(it)->pft != it->pft)
                    cleaned.push_back(*it);
            // delete entries that will never be used by branch and bound
            while (cleaned.size() > 1 && cleaned[cleaned.size() - 2].wt < min_lower[i])
                cleaned.pop_back();
            cleaned.push_back(LowerTableEntry{0, 0});
            knap.clear();
            for (auto it = std::prev(cleaned.end()); it != cleaned.begin(); --it)
                if (std::prev(it)->wt != it->wt)
                    knap.push_back(*it);
            knap.push_back(cleaned.front());
            knap.shrink_to_fit();
        }

        // merge adjacent equivalent columns
        std::vector<UpperTableEntry> cleaned;
        cleaned.reserve(cur.size());
        for (auto it = std::prev(cur.end()); it != cur.begin(); it--)
            if (!std::equal(std::prev(it)->L.begin(), std::prev(it)->L.end(), it->L.begin(), it->L.end()))
                cleaned.push_back(*it);
        // delete entries that will never be used by branch and bound
        while (cleaned.size() > 1 && cleaned[cleaned.size() - 2].wt < min_upper[i])
            cleaned.pop_back();
        cleaned.push_back(*cur.begin());
        cur.clear();
        std::reverse_copy(cleaned.begin(), cleaned.end(), std::back_inserter(cur));
        cur.shrink_to_fit();

        size_t n_entries = 0;
        for (auto col = cur.begin(); col != cur.end(); col++)
            n_entries += col->L.size();
        double density = n_entries/double(inst.up_cap+1)/double(inst.lo_cap+1);
        if (log) std::cerr << "sparse table for item " << i << " has n="
            << n_entries << ", density=" << density << std::endl;
        if (density >= density_threshold && i > 0) {
            if (log) std::cerr << "switching to dense table at item " << i << std::endl;
            compute_lower_bound_dense(i-1);
            break;
        }
    }
}


template<typename P>
void Comb_BKPSolver<P>::compute_lower_bound_dense(size_t start_at) {
    d_tab = (P*) malloc((start_at+2)*uint64_t(inst.up_cap+1)*(inst.lo_cap+1)*sizeof(P));
    ASSERT(d_tab && "couldn't allocate memory");
    dense_starts_at = start_at;

    omp_set_num_threads(num_threads);
    if (start_at == inst.n-1) {
        // if we are doing dense from the beginning
        // just initialize with zeros
        for (int ru = 0; ru <= inst.up_cap; ru++)
            for (int rl = 0; rl <= inst.lo_cap; rl++)
                d_tab_val(inst.n, ru, rl) = 0;
    } else {
        // initialize from the sparse table
        auto &tab = get_item_table_sparse(start_at+1);
        for (auto ue = tab.begin(); ue != tab.end(); ue++) {
            for (auto le = ue->L.begin(); le != ue->L.end(); le++) {
                if (std::next(le) != ue->L.end())
                    for (int j = le->wt; j < std::next(le)->wt; j++)
                        d_tab_val(start_at+1, ue->wt, j) = le->pft;
                else
                    for (int j = le->wt; j <= inst.lo_cap; j++)
                        d_tab_val(start_at+1, ue->wt, j) = le->pft;
            }
            if (std::next(ue) != tab.end())
                for (int k = ue->wt+1; k < std::next(ue)->wt; k++)
                    for (int j = 0; j <= inst.lo_cap; j++)
                        d_tab_val(start_at+1, k, j) = d_tab_val(start_at+1, ue->wt, j);
            else
                for (int k = ue->wt+1; k <= inst.up_cap; k++)
                    for (int j = 0; j <= inst.lo_cap; j++)
                        d_tab_val(start_at+1, k, j) = d_tab_val(start_at+1, ue->wt, j);
        }
    }
    #pragma omp parallel
    for (int i = start_at; i >= 0; i--) {
        #pragma omp barrier
        #pragma omp for nowait
        for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
            for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                d_tab_val(i, ru, rl) = std::min(d_tab_val(i+1, ru-inst.up_wt[i], rl),
                        std::max(P(d_tab_val(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i]),
                            d_tab_val(i+1, ru, rl)));

        if (min_lower[i] <= inst.lo_wt[i] && min_upper[i] <= inst.up_wt[i])
            #pragma omp for nowait
            for (int ru = 0; ru < inst.up_wt[i]; ru++)
                for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                    d_tab_val(i, ru, rl) = d_tab_val(i+1, ru, rl);

        if (min_lower[i] <= inst.lo_wt[i])
            #pragma omp for nowait
            for (int ru = std::max(min_upper[i], inst.up_wt[i]); ru <= inst.up_cap; ru++)
                for (int rl = 0; rl < inst.lo_wt[i]; rl++)
                    d_tab_val(i, ru, rl) = std::min(d_tab_val(i+1, ru-inst.up_wt[i], rl),
                            d_tab_val(i+1, ru, rl));

        if (min_upper[i] <= inst.up_wt[i])
            #pragma omp for nowait
            for (int ru = 0; ru < inst.up_wt[i]; ru++)
                for (int rl = std::max(min_lower[i], inst.lo_wt[i]); rl <= inst.lo_cap; rl++)
                    d_tab_val(i, ru, rl) = std::max(
                            P(d_tab_val(i+1, ru, rl-inst.lo_wt[i])+inst.pft[i]),
                            d_tab_val(i+1, ru, rl));
    }
}

template<typename P>
P &Comb_BKPSolver<P>::d_tab_val(size_t i, int up_cap, int lo_cap) {
    return d_tab[i*uint64_t(inst.up_cap+1)*(inst.lo_cap+1) + up_cap*uint64_t(inst.lo_cap+1) + lo_cap];
}

template<typename P>
std::vector<typename Comb_BKPSolver<P>::UpperTableEntry>& Comb_BKPSolver<P>::get_item_table_sparse(size_t i) {
    return s_tab[i];
}

template<typename P>
typename Comb_BKPSolver<P>::UpperTableEntry& Comb_BKPSolver<P>::get_upper_entry_sparse(size_t i, int up_cap) {
    auto &tab = get_item_table_sparse(i);
    auto uit = up_cap ? std::prev(std::lower_bound(tab.begin(), tab.end(), up_cap+1,
        [](const UpperTableEntry &element, int value){
            return element.wt < value;
    })) : tab.begin();
    ASSERT(uit->wt <= up_cap);
    return *uit;
}

template<typename P>
P Comb_BKPSolver<P>::s_tab_val(size_t i, int up_cap, int lo_cap) {
    UpperTableEntry& uit = get_upper_entry_sparse(i, up_cap);
    if (lo_cap == 0)
        return uit.L[0].pft;
    auto lit = std::prev(std::lower_bound(uit.L.begin(), uit.L.end(), lo_cap+1,
        [](const LowerTableEntry &element, int value){
            return element.wt < value;
    }));
    ASSERT(lit->wt <= lo_cap);
    return lit->pft;
}

template<typename P>
P Comb_BKPSolver<P>::lower_bound(size_t i, int up_cap, int lo_cap) {
    if (is_sparse[i])
        return s_tab_val(i, up_cap, lo_cap);
    else // dense
        return d_tab_val(i, up_cap, lo_cap);
}

template<typename P>
BKPSolution<max_pft_t> Comb_BKPSolver<P>::greedy_ub() {
    BKPSolution<max_pft_t> bkpsol(inst.n);

    knapsacksolver::Instance kpinst;
    kpinst.set_capacity(inst.up_cap);
    for (int i = 0; i < inst.n; i++)
        kpinst.add_item(inst.up_wt[i], inst.pft[i]);
    knapsacksolver::DynamicProgrammingPrimalDualOptionalParameters param;
    param.set_combo();
    auto output = knapsacksolver::dynamic_programming_primal_dual(kpinst, param);
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
    auto foutput = knapsacksolver::dynamic_programming_primal_dual(finst, param);
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

template<typename P>
void Comb_BKPSolver<P>::dump_lower_bound() {
    for (int i = 0; i < inst.n; i++) {
        std::ofstream os("lbdump_" + std::to_string(i) + ".pgm",
                       std::ios::out | std::ios::binary);
        os << "P2\n";
        os << (inst.up_cap+1) << " " << (inst.lo_cap+1) << "\n";
        P max_val = 1;
        for (int rl = 0; rl <= inst.lo_cap; rl++)
            for (int ru = 0; ru <= inst.up_cap; ru++)
                max_val = std::max(max_val, lower_bound(i, ru, rl));
        os << max_val << "\n";
        for (int rl = 0; rl <= inst.lo_cap; rl++) {
            for (int ru = 0; ru <= inst.up_cap; ru++)
                os << int(lower_bound(i, ru, rl)) << " ";
            os << "\n";
        }
        os.flush();
        os.close();
    }
}
