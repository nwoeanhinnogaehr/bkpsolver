# bkpsolver

This repo contains the code accompanying the paper titled "A Fast Combinatorial Algorithm for the Bilevel Knapsack Problem with Interdiction Constraints" by Noah Weninger and Ricardo Fukasawa.

Preprint: https://optimization-online.org/2022/10/a-fast-combinatorial-algorithm-for-the-bilevel-knapsack-problem-with-interdiction-constraints/

# Building

This code has only been tested on Linux. First, you will want to clone this repo.
The dependenceies are CMake, Gurobi, and [Knapsack Solver](https://github.com/fontanf/knapsacksolver).
CMake should be available from your distribution's package manager. Gurobi can be installed in any way you like, as long as `$GUROBI_HOME` is set appropriately. To install the knapsack solver, simply follow the build instructions from the [repo](https://github.com/fontanf/knapsacksolver), and then copy the built libraries into the `3rdparty/lib/` directory of this repo (overwriting the binaries that are there). If you run into linker issues, you may also need to copy the knapsacksolver header files into `3rdparty/include/knapsacksolver/`.

Now, create a build directory in the root of this repo:

``` $ mkdir build ```

Configure and build:

```
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

There is also a simple testing script included called `test.py`, which can be used to solve many instances and record timing information into a CSV file. To use the script, you can edit the code at the bottom of the file to use the desired solvers, tests, and CSV file.

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
