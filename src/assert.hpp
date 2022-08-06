#pragma once

#include <iostream>
#include <cstdlib>

#define ASSERT(x) if (!(x)) { \
    std::cerr << "ASSERTION FAILED: [[ " << #x << " ]] at " << __FILE__ << ":" << __LINE__;\
    exit(1); \
}
