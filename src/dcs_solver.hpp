#pragma once

#include "bkp_solver.hpp"

struct DCSTuple {
    int ptau = 0;
    int wtau = 0;
    std::vector<int> tau;
    bool added = false;
};

class DCS_BKPSolver : public BKPSolver<max_pft_t> {
public:
    DCS_BKPSolver(BKPInstance inst) : BKPSolver<max_pft_t>(inst), grb_env(true) { }

    BKPSolution<max_pft_t> solve() override;

    int alpha = 2, beta = 2, delta = 30, mu = 250, gamma = 5, omega = 3;
    int num_threads = 4;

private:
    double solve_model_obj(GRBModel &model);
    bool solve_model(GRBModel &model, BKPSolution<max_pft_t>& sol);
    void generate_tuples(int c);
    void add_tuple(DCSTuple &t, GRBModel &model, int c);
    void add_n_tuples(int c, int n, int weight_lim);
    GRBModel formulateCRIT2(size_t c, char type);
    bool NCR(BKPSolution<max_pft_t> &sol);
    void optimize_lower_level(BKPSolution<max_pft_t>& sol);
    void add_avoid_constraint(GRBModel &model, BKPSolution<max_pft_t> &sol);

    int cub, clb;
    GRBEnv grb_env;
    std::vector<GRBModel> crit2_lin, crit2_int;
    std::vector<std::vector<DCSTuple>> tuples;
    int ncons_added = 0;
};
