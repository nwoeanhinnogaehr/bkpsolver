# bkpsolver

This repo contains the code accompanying the paper titled "A Fast Combinatorial Algorithm for the Bilevel Knapsack Problem with Interdiction Constraints" by Noah Weninger and Ricardo Fukasawa.

Preprint: https://optimization-online.org/2022/10/a-fast-combinatorial-algorithm-for-the-bilevel-knapsack-problem-with-interdiction-constraints/

IPCO paper (conference version): https://doi.org/10.1007/978-3-031-32726-1_31

# Building

This code has only been tested on Linux. First, you will want to clone this repo.
The dependenceies are CMake, Gurobi, and [Knapsack Solver](https://github.com/fontanf/knapsacksolver).
CMake should be available from your distribution's package manager. Gurobi can be installed in any way you like, as long as `$GUROBI_HOME` is set appropriately. To install the knapsack solver, follow the build instructions from the [repo](https://github.com/fontanf/knapsacksolver), and then copy the following files libraries into the `3rdparty/lib/` directory of this repo (overwriting the binaries that are there):
- bazel-bin/knapsacksolver/libknapsacksolver.so
- bazel-bin/knapsacksolver/algorithms/libdynamic_programming_primal_dual.so
- bazel-bin/knapsacksolver/algorithms/libgreedy.so
- bazel-bin/knapsacksolver/algorithms/libsurrogate_relaxation.so
- bazel-bin/knapsacksolver/algorithms/libupper_bound_dantzig.so
- bazel-bin/external/optimizationtools/optimizationtools/utils/libinfo.so

The code was tested to work with commit 5464348be438e0b339f30c5f4f72cdaf701c99ec of knapsacksolver. Newer versions may work as well, but you may also need to copy the knapsacksolver header files into `3rdparty/include/knapsacksolver/` or fix the code accordingly for any changes made in knapsacksolver.

Now, create a build directory in the root of this repo:

``` $ mkdir build ```

Configure and build:

```
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
$ make
```

This should produce an executable in the `build` directory called `bkpsolver`.
Setting the compiler to clang is optional, but in my testing clang seems to do a better job.

# Usage

Instances must be in a specific file format. Many instances (with extension `.ki`) are
included in the `instances` directory. Details about the instance format are included in the next section.

Typical usage:

```
$ build/bkpsolver comb -p -l0 instances/CCLW/CCLW_n55_m2.ki # run Comb with standard parameters
$ build/bkpsolver dcs -a2 -b2 -g5 -d20 -m1000 -o5 instances/CCLW/CCLW_n55_m2.ki # run DCS with standard parameters
```

Full options:

```
SYNOPSIS
        ./bkpsolver dcs [-a <alpha=2>] [-b <beta=2>] [-g <gamma=5>] [-d <delta=30>] [-m <mu=250>]
                    [-o <omega=3>] [-q] [-j <num_threads=4>] [<file>]

        ./bkpsolver comb [-l <lookback=1>] [-p] [-w] [-lb-only] [-q] [-j <num_threads=4>] [<file>]

OPTIONS
        dcs         use algorithm by Federico Della Croce and Rosario Scatamacchia
        comb        use combinatorial algorithm
        -p          search for best prefix (try this if number of leaves is large)
        -q          quiet mode: do not log to stderr
        <file>      filename (will read from stdin if absent)
```

There is also a simple testing script included called `run_tests.py`, which can be used to solve many
instances and record timing information into a JSON file. To use the script, you can edit the code at
the bottom of the file to use the desired solvers, tests, and JSON file.

## Troubleshooting Gurobi
Bkpsolver currently targets Gurobi 10.0 and is not guaranteed to work with any other version. If you
recieve the error message
```terminate called after throwing an instance of 'GRBException'```
then your Gurobi licence file may be invalid or expired. Note that with an academic license you need
to reactivate Gurobi every year or whenever you upgrade to a new major version.

# Instance format
The instances are simple text files with extension `.ki`, in the following format. The file contains at least 6 lines.

1. The first line contains $n$, the number of items.
2. The second line contains $C^L$, the lower-level capacity.
3. The third line contains $C^U$, the upper-level capacity.
4. The fourth line contains $n$ integers $w^L_1,\dots,w^L_n$ representing the lower-level weights.
5. The fifth line contains $n$ integers $w^U_1,\dots,w^U_n$ representing the upper-level weights.
6. The sixth line contains $n$ integers $p_1,\dots,p_n$ representing the profits.

After these 6 lines, the remaining lines of the file may contain metadata information. This information is ignored by the solver but is used by the testing script.

The accompanying answer file (for instances where the solution is known) has extension `.ans` and contains a single number: the profit value of an optimal solution.

The file `instances_lit.tar.zst` contains a collection of instances from other papers in the literature, in the original format they were released in. It also contains scripts to convert each of these formats into the simple text format that we use.

# Output format

Here is an example of some output from the solver, annotated with the meaning of each line.
All time measurements are wall-clock time, in milliseconds.

```
algorithm comb # which algorithm was used (either comb or dcs)
file instances/CCLW/CCLW_n35_m0.ki # filename of instance
threads 4 # number of threads used in computation
lookback 0 # how far to look back when greedily computing LB prefix (not used in paper)
best_prefix 1 # whether to use the best prefix (always set to 1 in paper)
weak_lb 0 # whether to use the weak lower bound
lb_only 0 # whether to only compute the lower bound
greedy_ub 285 # profit obtained by GreedyHeuristic upper bound
dcs_lb 218 # LP based lower bound, similar to DCS lower bound
initial_bound_test_time 11.5848 # time in milliseconds to compute greedy_ub and dcs_lb
dp_lb 279 # DP lower bound value
dp_lb_time 0.915699 # time to compute DP lower bound
nodes 45 # number of branch and bound nodes explored
leaves 1 # number of branch and bound leaves
bnb_time 0.28387 # total time for branch and bound
opt_time 0.272626 # time for branch and bound to find an optimal solution
proof_time 0.011244 # time for branch and bound to prove optimality
total_time 12.8214 # total time to solve the instance
max_depth 35 # max depth explored by branch and bound
avg_depth 15.1556 # average depth explored by branch and bound
pft_bits 16 # number of bits used to store profit values
memory 23.9492 # maximum memory usage in megabytes
upper 00000000100100001001100000001001000 # upper level solution (item 1 on left)
lower 10101000000000100000010000000000100 # lower level solution 
up_util 152 # how much upper level capacity was used
lo_util 155 # how much lower level capacity was used
profit 279 # solution profit value
```
