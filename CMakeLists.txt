cmake_minimum_required(VERSION 3.3)

project(BKPSolver)

add_compile_options(
    -march=native
    -mtune=native
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(GUROBI REQUIRED)
find_package(OpenMP)

include_directories(${CMAKE_SOURCE_DIR}/3rdparty/include)
include_directories(${GUROBI_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/3rdparty/lib)

add_executable(bkpsolver
   src/bkp_instance.cpp
   src/dcs_solver.cpp
   src/main.cpp)
target_link_libraries(bkpsolver knapsacksolver minknap surrelax greedy dantzig info)
target_link_libraries(bkpsolver gurobi_g++5.2)
target_link_libraries(bkpsolver ${GUROBI_CXX_LIBRARY})
target_link_libraries(bkpsolver ${GUROBI_LIBRARY})
target_link_libraries(bkpsolver OpenMP::OpenMP_CXX)
