#pragma once

#include <vector>
#include <cstddef>
#include "assert.hpp"

class BKPInstance {
public:
    size_t n;
    std::vector<int> lo_wt, up_wt, pft;
    int lo_cap, up_cap;

    BKPInstance() { }

    BKPInstance(size_t n, std::vector<int> pft, std::vector<int> lo_wt,
                std::vector<int> up_wt, int lo_cap, int up_cap)
        : n(n), pft(pft), lo_wt(lo_wt), up_wt(up_wt), lo_cap(lo_cap), up_cap(up_cap) {
        ASSERT(n >= 1);
        ASSERT(pft.size() == n);
        ASSERT(lo_wt.size() == n);
        ASSERT(up_wt.size() == n);
        ASSERT(lo_cap >= 0);
        ASSERT(up_cap >= 0);
        for(size_t i = 0; i < n; i++) {
            ASSERT(pft[i] >= 0);
            ASSERT(lo_wt[i] >= 0);
            ASSERT(up_wt[i] >= 0);
        }
    }

    static BKPInstance from_stdin();
    static BKPInstance from_file(std::string filename);
};
