#include "bkp_instance.hpp"
#include "bkp_solver.hpp"
#include "dcs_solver.hpp"

#include "knapsacksolver/solution.hpp"
#include "knapsacksolver/algorithms/minknap.hpp"
#include <gurobi_c++.h>

#include <bits/stdc++.h>
using namespace std;

BKPSolution<max_pft_t> DCS_BKPSolver::solve() {
    grb_env.set("OutputFlag", "0");
    grb_env.set("MIPFocus", "2");
    grb_env.set("Presolve", "2");
    grb_env.set("PreSparsify", "1");
    grb_env.set("Cuts", "0");
    grb_env.set("Threads", to_string(num_threads));
    grb_env.set("MIPGap", "0");
    grb_env.start();

    sort_inst();
    if (log) cerr << "computing NCR..." << endl;
    BKPSolution<max_pft_t> ub(inst.n);
    NCR(ub);
    if (log) cerr << "computing initial lower bounds..." << endl;
    clb = find_crit_lb();
    cub = min(int(inst.n)-1, find_crit_ub());
    tuples.resize(inst.n);
    vector<double> lp_vals;
    for (int c = clb; c <= cub; c++) {
        generate_tuples(c);
        crit2_lin.push_back(formulateCRIT2(c, GRB_CONTINUOUS));
        crit2_int.push_back(formulateCRIT2(c, GRB_BINARY));
        add_n_tuples(c, mu, inst.lo_cap);
        double val = solve_model_obj(crit2_lin.back());
        lp_vals.push_back(val);
    }
    vector<int> L(cub-clb+1);
    iota(L.begin(), L.end(), 0);
    sort(L.begin(), L.end(), [&](int a, int b) {
        return lp_vals[a] < lp_vals[b];
    });
    if (log && L.size()) cerr << "lower bound: " << lp_vals[L[0]] << endl;
    if (L.size() == 0 || lp_vals[L[0]] >= ub.pft) {
        unsort_sol(ub);
        return ub;
    }
    if (log) cerr << "computing initial upper bound..." << endl;
    for (int i = 0; i < min(int(L.size()), gamma); i++)
        if (lp_vals[L[i]] < ub.pft) {
            BKPSolution<max_pft_t> sol(inst.n);
            if (solve_model(crit2_int[L[i]], sol))
                if (sol.pft < ub.pft) {
                    optimize_lower_level(sol);
                    if (sol.pft < ub.pft)
                        ub = sol;
                }
        }
    if (log) cerr << "upper bound: " << ub.pft << endl;
    for (int i = 0; i < L.size(); i++) {
        if (log) cerr << "step " << i << " (c=" << L[i]+clb << ")" << endl;
        if (lp_vals[L[i]] >= ub.pft) {
            unsort_sol(ub);
            return ub;
        }
        solve_model_obj(crit2_lin[L[i]]);
        int nfixed = 0;
        for (int j = 0; j < inst.n; j++) {
            GRBVar xj = crit2_lin[L[i]].getVarByName("x" + to_string(j));
            double rc = xj.get(GRB_DoubleAttr_RC);
            if (abs(rc) >= ub.pft - lp_vals[L[i]]) {
                GRBVar xj_int = crit2_int[L[i]].getVarByName("x" + to_string(j));
                crit2_int[L[i]].addConstr(xj_int == xj.get(GRB_DoubleAttr_X), "fix" + to_string(ncons_added++));
                nfixed++;
            }
        }
        for (int j = 0; j < inst.lo_wt[L[i]+clb]; j++) {
            GRBVar kj = crit2_lin[L[i]].getVarByName("k" + to_string(j));
            double rc = kj.get(GRB_DoubleAttr_RC);
            if (abs(rc) >= ub.pft - lp_vals[L[i]]) {
                GRBVar kj_int = crit2_int[L[i]].getVarByName("k" + to_string(j));
                crit2_int[L[i]].addConstr(kj_int == kj.get(GRB_DoubleAttr_X), "fix" + to_string(ncons_added++));
                nfixed++;
            }
        }
        if (nfixed)
            if (log) cerr << "fixed " << nfixed << " vars" << endl;
        crit2_int[L[i]].update();
        crit2_int[L[i]].set("Cutoff", to_string(ub.pft));
        BKPSolution<max_pft_t> sol(inst.n);
        if (!solve_model(crit2_int[L[i]], sol))
            continue;
        int bestlb = sol.pft;
        int niters_same = 1;
        while (sol.pft < ub.pft) {
            int lb = sol.pft;
            if (log) cerr << "lower bound: " << lb << endl;
            optimize_lower_level(sol);
            if (sol.pft < ub.pft) {
                ub = sol;
                if (log) cerr << "new upper bound: " << ub.pft << endl;
                if (lb >= ub.pft)
                    break;
            }
            DCSTuple tup;
            for (int j = 0; j < L[i]+clb; j++)
                if (!sol.up_sol[j] && !sol.lo_sol[j]) {
                    tup.tau.push_back(j);
                    tup.ptau -= inst.pft[j];
                    tup.wtau -= inst.lo_wt[j];
                }
            for (int j = L[i]+clb; j < inst.n; j++)
                if (!sol.up_sol[j] && sol.lo_sol[j]) {
                    tup.tau.push_back(j);
                    tup.ptau += inst.pft[j];
                    tup.wtau += inst.lo_wt[j];
                }
            if (tup.tau.size() && tup.ptau > 0 && tup.wtau <= inst.lo_wt[clb+L[i]]) {
                if (log) {
                    cerr << "add tuple cons: ";
                    for (int i : tup.tau)
                        cerr << i << " ";
                    cerr << endl;
                }
                add_tuple(tup, crit2_int[L[i]], clb+L[i]);
            } else {
                if (log) {
                    cerr << "add avoid cons: ";
                    for (size_t i = 0; i < inst.n; i++)
                        cerr << (int)sol.lo_sol[i];
                    cerr << endl;
                }
                add_avoid_constraint(crit2_int[L[i]], sol);
            }
            if (niters_same == omega) {
                niters_same = 0;
                int weight_lim = inst.lo_cap;
                for (int i = 0; i < L[i]+clb; i++)
                    weight_lim -= inst.lo_wt[i]*(1 - sol.up_sol[i]);
                add_n_tuples(L[i]+clb, mu, weight_lim);
            }
            crit2_int[L[i]].update();
            crit2_int[L[i]].set("Cutoff", to_string(ub.pft));
            solve_model(crit2_int[L[i]], sol);
            if (sol.pft <= bestlb)
                niters_same++;
            else
                niters_same = 1;
            bestlb = max(bestlb, (int)sol.pft);
        }
    }

    unsort_sol(ub);
    return ub;
}
void DCS_BKPSolver::add_avoid_constraint(GRBModel &model, BKPSolution<max_pft_t> &sol) {
    GRBLinExpr expr = 0.0;
    for (int i = 0; i < inst.n; i++)
        if (sol.lo_sol[i])
            expr += model.getVarByName("x" + to_string(i));
    model.addConstr(expr >= 1, "avoid" + to_string(ncons_added++));
}
double DCS_BKPSolver::solve_model_obj(GRBModel &model) {
    model.optimize();
    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
        return INFINITY;
    return model.get(GRB_DoubleAttr_ObjVal);
}
bool DCS_BKPSolver::solve_model(GRBModel &model, BKPSolution<max_pft_t>& sol) {
    model.optimize();
    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
        return false;
    ASSERT(int(0.5+model.get(GRB_DoubleAttr_ObjVal)) < (1ull<<(8*sizeof(max_pft_t))));
    sol.pft = int(0.5+model.get(GRB_DoubleAttr_ObjVal));
    for (size_t i = 0; i < inst.n; i++) {
        sol.up_sol[i] = model.getVarByName("x" + to_string(i)).get(GRB_DoubleAttr_X);
        sol.lo_sol[i] = 1-sol.up_sol[i];
    }
    return true;
}
void DCS_BKPSolver::optimize_lower_level(BKPSolution<max_pft_t>& sol) {
    knapsacksolver::Instance finst;
    finst.set_capacity(inst.lo_cap);
    for (int i = 0; i < inst.n; i++) {
        if (!sol.up_sol[i])
            finst.add_item(inst.lo_wt[i], inst.pft[i]);
        else
            finst.add_item(inst.lo_wt[i], 0);
    }
    knapsacksolver::MinknapOptionalParameters param;
    param.set_combo();
    auto foutput = knapsacksolver::minknap(finst, param);
    auto fsol = foutput.solution;
    for (int i = 0; i < inst.n; i++) {
        if (sol.up_sol[i]) sol.lo_sol[i] = 0;
        else if (fsol.contains_idx(i)) sol.lo_sol[i] = 1;
        else sol.lo_sol[i] = 0;
    }
    sol.pft = fsol.profit();
}
void DCS_BKPSolver::generate_tuples(int c) {
    int a = max(c-delta, 0);
    int b = min(c+delta, int(inst.n)-1);
    vector<vector<int>> backward_sets, forward_sets;
    backward_sets.emplace_back();
    for (int i = a; i < c; i++) {
        size_t sz = backward_sets.size();
        for (size_t j = 0; j < sz; j++)
            if (backward_sets[j].size() < alpha) {
                vector<int> s = backward_sets[j];
                s.push_back(i);
                backward_sets.push_back(s);
            }
    }
    forward_sets.emplace_back();
    for (int i = c; i <= b; i++) {
        size_t sz = forward_sets.size();
        for (size_t j = 0; j < sz; j++)
            if (forward_sets[j].size() < beta) {
                vector<int> s = forward_sets[j];
                s.push_back(i);
                forward_sets.push_back(s);
            }
    }
    for (vector<int> &f : forward_sets)
        for (vector<int> &b : backward_sets) {
            int ptau = 0, wtau = 0;
            for (int i : b) ptau -= inst.pft[i];
            for (int i : f) ptau += inst.pft[i];
            if (ptau <= 0) continue;
            for (int i : b) wtau -= inst.lo_wt[i];
            for (int i : f) wtau += inst.lo_wt[i];
            if (wtau > inst.lo_wt[c]) continue;
            tuples[c].emplace_back();
            tuples[c].back().ptau = ptau;
            tuples[c].back().wtau = wtau;
            for (int i : b) tuples[c].back().tau.push_back(i);
            for (int i : f) tuples[c].back().tau.push_back(i);
        }
    sort(tuples[c].begin(), tuples[c].end(), [](const DCSTuple &a, const DCSTuple &b) {
        return a.tau.size() == b.tau.size() ?
            a.ptau*b.wtau > b.ptau*a.wtau : a.tau.size() < b.tau.size();
    });
}
void DCS_BKPSolver::add_tuple(DCSTuple &t, GRBModel &model, int c) {
    GRBVar var_pi = model.getVarByName("pi");
    GRBLinExpr expr = 0.0;
    for (int i : t.tau)
        expr -= model.getVarByName("x" + to_string(i));
    if (t.wtau <= 0)
        expr += 1;
    else
        for (int i = t.wtau; i <= inst.lo_wt[c]; i++)
            expr += model.getVarByName("k" + to_string(i-1));
    expr *= t.ptau;
    model.addConstr(var_pi >= expr, "tup" + to_string(ncons_added++));
}
void DCS_BKPSolver::add_n_tuples(int c, int n, int weight_lim) {
    int nadded = 0;
    for (int i = 0; i < tuples[c].size() && nadded < n; i++) {
        if (tuples[c][i].added || tuples[c][i].wtau > weight_lim) continue;
        add_tuple(tuples[c][i], crit2_lin[c-clb], c);
        add_tuple(tuples[c][i], crit2_int[c-clb], c);
        nadded++;
        tuples[c][i].added = true;
    }
    crit2_int[c-clb].update();
    crit2_lin[c-clb].update();
    if (nadded)
        if (log) cerr << "added " << nadded << " tuples for c=" << c << endl;
}
GRBModel DCS_BKPSolver::formulateCRIT2(size_t c, char type) {
    GRBModel model(grb_env);

    vector<GRBVar> vars_x, vars_k;
    GRBVar var_pi = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "pi");
    GRBLinExpr obj = var_pi;
    for (size_t i = 0; i < inst.n; i++) {
        GRBVar var = model.addVar(0.0, 1.0, 0.0, type, "x" + to_string(i));
        vars_x.push_back(var);
        if (i < c)
            obj += inst.pft[i] - var*inst.pft[i];
    }
    for (int i = 0; i < inst.lo_wt[c]; i++)
        vars_k.push_back(model.addVar(0.0, 1.0, 0.0, type, "k" + to_string(i)));
    model.setObjective(obj, GRB_MINIMIZE);

    GRBLinExpr up_kp = 0.0;
    for (size_t i = 0; i < inst.n; i++)
        up_kp += vars_x[i]*inst.up_wt[i];
    model.addConstr(up_kp <= inst.up_cap, "up_kp");
    GRBLinExpr lo_kp = 0.0;
    for (size_t i = 0; i < c; i++)
        lo_kp += inst.lo_wt[i] - vars_x[i]*inst.lo_wt[i];
    for (int i = 0; i < inst.lo_wt[c]; i++)
        lo_kp += (i+1)*vars_k[i];
    model.addConstr(lo_kp == inst.lo_cap, "lo_kp");
    GRBLinExpr sum_k = 0.0;
    for (int i = 0; i < inst.lo_wt[c]; i++)
        sum_k += vars_k[i];
    model.addConstr(sum_k == 1, "one_K");
    model.addConstr(vars_x[c] == 0, "xc_zero");
    model.addConstr(var_pi >= inst.pft[c]*vars_k[inst.lo_wt[c]-1], "try_crit");

    model.update();

    return model;
}
bool DCS_BKPSolver::NCR(BKPSolution<max_pft_t> &sol) {
    GRBModel model(grb_env);

    vector<GRBVar> vars;
    GRBLinExpr obj = 0.0;
    for (size_t i = 0; i < inst.n; i++) {
        GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
        vars.push_back(var);
        obj += inst.pft[i] - var*inst.pft[i];
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBLinExpr up_kp = 0.0;
    for (size_t i = 0; i < inst.n; i++)
        up_kp += vars[i]*inst.up_wt[i];
    model.addConstr(up_kp <= inst.up_cap, "up_kp");
    GRBLinExpr lo_kp = 0.0;
    for (size_t i = 0; i < inst.n; i++)
        lo_kp += inst.lo_wt[i] - vars[i]*inst.lo_wt[i];
    model.addConstr(lo_kp <= inst.lo_cap, "lo_kp");

    model.optimize();
    if (model.get(GRB_IntAttr_SolCount) == 0) {
        sol.pft = (max_pft_t)-1;
        return false;
    }
    ASSERT(int(0.5+model.get(GRB_DoubleAttr_ObjVal)) < (1ull<<(8*sizeof(max_pft_t))));
    sol.pft = int(0.5+model.get(GRB_DoubleAttr_ObjVal));
    for (size_t i = 0; i < inst.n; i++) {
        sol.up_sol[i] = vars[i].get(GRB_DoubleAttr_X);
        sol.lo_sol[i] = 1-vars[i].get(GRB_DoubleAttr_X);
    }
    return true;
}
