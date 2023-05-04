#include "bkp_instance.hpp"
#include <iostream>
#include <fstream>

using namespace std;

BKPInstance BKPInstance::from_stdin() {
    size_t n;
    vector<int> lo_wt, up_wt, pft, cnt;
    int lo_cap, up_cap;

    cin >> n;
    pft.resize(n);
    lo_wt.resize(n);
    up_wt.resize(n);
    cnt.resize(n, 1);
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
    string s;
    std::getline(cin, s);
    if (cin.peek() != EOF)
        cin >> s;
    if (s == "cnt")
        for(int i = 0; i < n; i++)
            cin >> cnt[i];
    return BKPInstance(n, pft, lo_wt, up_wt, cnt, lo_cap, up_cap);
}

BKPInstance BKPInstance::from_file(std::string filename) {
    size_t n;
    vector<int> lo_wt, up_wt, pft, cnt;
    int lo_cap, up_cap;

    ifstream f(filename);
    ASSERT(f.good());
    f >> n;
    pft.resize(n);
    lo_wt.resize(n);
    up_wt.resize(n);
    cnt.resize(n, 1);
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
    string s;
    std::getline(f, s);
    if (f.peek() != EOF)
        f >> s;
    if (s == "cnt")
        for(int i = 0; i < n; i++)
            f >> cnt[i];
    return BKPInstance(n, pft, lo_wt, up_wt, cnt, lo_cap, up_cap);
}
