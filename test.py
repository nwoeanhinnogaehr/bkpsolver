import os
import subprocess
import time

class Result:
    def __init__(self):
        self.error = None
        self.err_extra = ""
        self.opt = None
        self.cputime = 0
        self.expected_opt = 0

    def is_success(self):
        return self.error is None

    def __str__(self):
        out = ""
        if self.is_success():
            out = "SUCCESS"
        else:
            out = "ERROR({}): {}".format(self.error, self.err_extra)
        out = out + " Opt={} Expect={} Time={:.3f}".format(self.opt, self.expected_opt, self.cputime)
        return out

GENERATE_SOL = "GENERATE_SOL"
EXACT = "EXACT"
LOWER_BOUND = "LOWER_BOUND"
UPPER_BOUND = "UPPER_BOUND"

class Solver:
    def __init__(self, filename, args, mode, approx=None):
        self.filename = filename
        self.args = args
        self.mode = mode
        assert(approx is None) # TODO
        self.approx = approx

    def __str__(self):
        return "\"{} {}\" ({})".format(self.filename, " ".join(self.args), self.mode)

ERR_TIMEOUT = "timeout"
ERR_MEMORY = "memory"
ERR_RETCODE = "retcode"
ERR_PARSE = "parse"
ERR_WRONG = "wrong answer"
ERR_OTHER = "other"

def run_test(solver, test, timeout):
    print("Testing {} on {}".format(solver, test))
    print("    ", end="")
    result = Result()

    if solver.mode == GENERATE_SOL:
        if os.path.exists(test[:-2] + "ans"):
            print("answer file exists, not generating")
            return result

    with open(test, "r") as f:
        testdata = f.read()

    start_time = time.time()
    proc = subprocess.Popen([solver.filename, *solver.args],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            #stderr=subprocess.PIPE,
            universal_newlines=True)
    try:
        stdout_data, stderr_data = proc.communicate(input=testdata, timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        result.error = ERR_TIMEOUT
        result.cputime = timeout
        return result
    proc.kill()
    end_time = time.time()
    result.cputime = end_time - start_time

    if proc.returncode != 0:
        result.error = ERR_RETCODE
        result.err_extra = "{}".format(proc.returncode)
        return result

    try:
        result.opt = int(stdout_data.strip())
    except ValueError:
        result.error = ERR_PARSE
        result.err_extra = stdout_data.strip()
        return result

    try:
        with open(test[:-2] + "ans", "r") as f:
            result.expected_opt = int(f.read())
    except FileNotFoundError:
        result.expected_opt = result.opt
        print("answer file doesn't exist, generating")
        with open(test[:-2] + "ans", "w") as f:
            f.write(str(result.opt))
            f.write("\n")
    if solver.mode == EXACT:
        if result.expected_opt != result.opt:
            result.error = ERR_WRONG
            result.err_extra = "expected == {}".format(result.expected_opt)
            return result
    if solver.mode == LOWER_BOUND:
        if result.expected_opt < result.opt:
            result.error = ERR_WRONG
            result.err_extra = "expected <= {}".format(result.expected_opt)
            return result
    if solver.mode == UPPER_BOUND:
        if result.expected_opt > result.opt:
            result.error = ERR_WRONG
            result.err_extra = "expected >= {}".format(result.expected_opt)
            return result

    return result

def run_tests(solvers, tests, f, timeout=3):
    f.write("solver,instance,n,lo_cap,up_cap,group,class,ins,err,opt,cputime\n")
    nfail = 0
    for solver in solvers:
        for test in tests:
            result = run_test(solver, test.filename, timeout)
            print()
            print(result)
            if not result.is_success():
                nfail += 1
            f.write('{} {},{},{},{},{},{},{},{},{},{},{}\n'.format(
                    solver.filename,
                    " ".join(solver.args),
                    test.filename,
                    test.n,
                    test.lo_cap,
                    test.up_cap,
                    test.group,
                    test.cls,
                    test.ins,
                    result.error,
                    result.opt,
                    result.cputime
                    ))
            f.flush()
    if nfail > 0:
        print("\n{} FAILURES".format(nfail))

class Test:
    def __init__(self, filename):
        self.filename = filename
        with open(filename, "r") as f:
            testdata = f.read()
        lines = testdata.split("\n")
        self.n = int(lines[0].strip())
        self.lo_cap = int(lines[1].strip())
        self.up_cap = int(lines[2].strip())
        self.group = None
        self.cls = None
        self.ins = None
        for l in lines[6:]:
            if len(l.strip()) == 0:
                continue
            key, val = l.strip().split()
            if key == "group":
                self.group = val
            if key == "class":
                self.cls = val
            if key == "ins":
                self.ins = val

    def __str__(self):
        return "{} {} {} {} {} {} {}\n".format(self.filename, self.n, self.lo_cap, self.up_cap, self.group, self.cls, self.ins)

def get_all_tests():
    tests = []
    for group in sorted(os.listdir("instances")):
        files = os.listdir("instances/" + group)
        for file in sorted(files):
            if file[-3:] == ".ki":
                tests.append(Test("instances/" + group + "/" + file))
    return tests

def get_all_solvers():
    return [
            Solver("build/bkpsolver", ["comb"], EXACT),
            Solver("build/bkpsolver", ["comb", "-l", "0"], EXACT),
            Solver("build/bkpsolver", ["comb", "-p"], EXACT),
            Solver("build/bkpsolver", ["dcs"], EXACT),
            Solver("build/bkpsolver", ["dcs", "-a", "2", "-b", "2", "-g", "5", "-d", "20", "-m", "1000", "-o", "5"], EXACT)
    ]

tests = get_all_tests()
#tests = list(filter(lambda t: t.group == "CCLW", tests))
with open("results.csv", "w") as f:
    run_tests(get_all_solvers(), tests, f, timeout=900)
