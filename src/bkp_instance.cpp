#include "bkp_instance.hpp"
#include <iostream>
#include <fstream>

using namespace std;

BKPInstance BKPInstance::from_stdin() {
    size_t n;
    vector<int> lo_wt, up_wt, pft;
    int lo_cap, up_cap;

    cin >> n;
    pft.resize(n);
    lo_wt.resize(n);
    up_wt.resize(n);
    cin >> lo_cap;
    cin >> up_cap;
    for(int i = 0; i < n; i++) {
        cin >> lo_wt[i];
        ASSERT(lo_wt[i] <= lo_cap);
    }
    for(int i = 0; i < n; i++) {
        cin >> up_wt[i];
        ASSERT(up_wt[i] <= up_cap);
    }
    for(int i = 0; i < n; i++)
        cin >> pft[i];
    return BKPInstance(n, pft, lo_wt, up_wt, lo_cap, up_cap);
}

BKPInstance BKPInstance::from_file(std::string filename) {
    size_t n;
    vector<int> lo_wt, up_wt, pft;
    int lo_cap, up_cap;

    ifstream f(filename);
    ASSERT(f.good());
    f >> n;
    pft.resize(n);
    lo_wt.resize(n);
    up_wt.resize(n);
    f >> lo_cap;
    f >> up_cap;
    for(int i = 0; i < n; i++) {
        f >> lo_wt[i];
        ASSERT(lo_wt[i] <= lo_cap);
    }
    for(int i = 0; i < n; i++) {
        f >> up_wt[i];
        ASSERT(up_wt[i] <= up_cap);
    }
    for(int i = 0; i < n; i++)
        f >> pft[i];
    return BKPInstance(n, pft, lo_wt, up_wt, lo_cap, up_cap);
}
