import os
import subprocess
import time

GENERATE_SOL = "GENERATE_SOL"
EXACT = "EXACT"
LOWER_BOUND = "LOWER_BOUND"
UPPER_BOUND = "UPPER_BOUND"

class Solver:
    def __init__(self, filename, args, mode):
        self.filename = filename
        self.args = args
        self.mode = mode
    def __str__(self):
        return "\"{} {}\" ({})".format(self.filename, " ".join(self.args), self.mode)


class Result:
    def __init__(self):
        self.error = None
        self.err_extra = ""
        self.opt = None
        self.expected_opt = 0

    def is_success(self):
        return self.error is None

    def __str__(self):
        out = ""
        if self.is_success():
            out = "SUCCESS"
        else:
            out = "ERROR({}): {}".format(self.error, self.err_extra)
        out = out + " Opt={} Expect={}".format(self.opt, self.expected_opt)
        return out
 
ERR_TIMEOUT = "timeout"
ERR_MEMORY = "memory"
ERR_RETCODE = "retcode"
ERR_PARSE = "parse"
ERR_WRONG = "wronganswer"
ERR_OTHER = "other"

def run_test(solver, test, timeout):
    print("Testing {} on {}".format(solver, test.filename))
    print("    ", end="")
    result = {
        "n": str(test.n),
        "lo_cap": str(test.lo_cap),
        "up_cap": str(test.up_cap),
        "group": str(test.group),
        "cls": str(test.cls),
        "ins": str(test.ins),
        "status": "ok"
    }

    result["solver_cmdline"] = "_".join(solver.args)
    output_file = output_dir + "/" + test.filename[10:-2] + "_".join(solver.args) + ".result"

    if solver.mode == GENERATE_SOL:
        if os.path.exists(test.filename[:-2] + "ans"):
            print("answer file exists, not generating")
            return result
    else:
        if os.path.exists(output_file):
            result["status"] = "skip"
            return result
            
    proc = subprocess.Popen([solver.filename, *solver.args, test.filename],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            #stderr=subprocess.PIPE,
            universal_newlines=True)
    try:
        stdout_data, stderr_data = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout_data, stderr_data = proc.communicate()
        result["status"] = ERR_TIMEOUT
        result["total_time"] = timeout * 1000
    proc.kill()

    if result["status"] != ERR_TIMEOUT and proc.returncode != 0:
        result["status"] = ERR_RETCODE
        result["error"] = "{}".format(proc.returncode)

    try:
        for line in stdout_data.splitlines():
            key, val = line.strip().split()
            result[key] = val
    except ValueError:
        result["status"] = ERR_PARSE
        result["error"] = stdout_data.strip()
        return result

    if result["status"] == "ok":
        try:
            with open(test.filename[:-2] + "ans", "r") as f:
                result["expected_profit"] = f.read()
        except FileNotFoundError:
            print("answer file doesn't exist, generating")
            result["expected_profit"] = result["profit"].strip()
            with open(test.filename[:-2] + "ans", "w") as f:
                f.write(result["profit"])
                f.write("\n")

        if solver.mode == EXACT:
            if int(result["expected_profit"]) != int(result["profit"]):
                result["status"] = ERR_WRONG
                result["error"] = "expected=={}".format(result["expected_profit"])
        if solver.mode == LOWER_BOUND:
            if int(result["expected_profit"]) < int(result["profit"]):
                result["status"] = ERR_WRONG
                result["error"] = "expected<={}".format(result["expected_profit"])
        if solver.mode == UPPER_BOUND:
            if int(result["expected_profit"]) > int(result["profit"]):
                result["status"] = ERR_WRONG
                result["error"] = "expected>={}".format(result["expected_profit"])

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        for key in result:
            f.write(str(key) + " " + str(result[key]) + "\n")
    
    return result

def run_tests(solvers, tests, f, timeout=3):
    nfail = 0
    for solver in solvers:
        for test in tests:
            result = run_test(solver, test, timeout)
            print(result["status"])
            if result["status"] not in ["ok", "skip"]:
                nfail += 1
            if result["status"] != "skip":
                f.write(str(result) + "\n")
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

output_dir = "test_results"
all_solvers = [
    Solver("build/bkpsolver", ["comb", "-l0", "-p", "-q", "-j1"], EXACT),
    Solver("build/bkpsolver", ["comb", "-l0", "-p", "-q", "-j4"], EXACT),
    Solver("build/bkpsolver", ["comb", "-l0", "-p", "-q", "-j16"], EXACT),
    #Solver("build/bkpsolver", ["dcs", "-a2", "-b2", "-g5", "-d20", "-m1000", "-o5", "-q"], EXACT)
]

tests = get_all_tests()
new_tests = list(filter(lambda t: t.group in ["New", None], tests))
lit_tests = list(filter(lambda t: t.group in ["CCLW", "DCS", "FMS", "TRS", "DeNegre"], tests))
os.makedirs(output_dir, exist_ok=True)
with open(output_dir + "/all_results", "a") as f:
    run_tests(all_solvers, lit_tests, f, timeout=3600)
    run_tests(all_solvers, new_tests, f, timeout=900)
