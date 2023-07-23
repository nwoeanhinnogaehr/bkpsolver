from random import randint, shuffle
import math

def uniform_array(n, lo, hi):
    return [randint(lo, hi) for i in range(n)]

def uniform_array_nodup(n, lo, hi, skip=1):
    assert(n <= (hi-lo+skip)//skip)
    a = list(range(lo, hi+1, skip))
    shuffle(a)
    return a[:n]

def add(x, y):
    return [z+y for z in x]

num_inst = 0
def write_instance(n, p, cu, cl, wu, wl, cnt, meta={"group":"Count"}):
    global num_inst
    with open("cnt{}.ki".format(num_inst), "w") as f:
        f.write("{}\n{}\n{}\n".format(n,cl,cu))
        f.write("{}\n".format(" ".join(map(str, wl))))
        f.write("{}\n".format(" ".join(map(str, wu))))
        f.write("{}\n".format(" ".join(map(str, p))))
        f.write("cnt {}\n".format(" ".join(map(str, cnt))))
        for key in meta:
            f.write("{} {}\n".format(key, meta[key]))
    with open("cln{}.ki".format(num_inst), "w") as f:
        cn = 0
        cwl = []
        cwu = []
        cp = []
        for i in range(n):
            cn += cnt[i]
            for j in range(cnt[i]):
                cwl.append(wl[i])
                cwu.append(wu[i])
                cp.append(p[i])
        f.write("{}\n{}\n{}\n".format(cn,cl,cu))
        f.write("{}\n".format(" ".join(map(str, cwl))))
        f.write("{}\n".format(" ".join(map(str, cwu))))
        f.write("{}\n".format(" ".join(map(str, cp))))
        for key in meta:
            f.write("{} {}\n".format(key, meta[key]))
        num_inst += 1

sizes = [5, 10, 15]
n_of_each = 1
R = [10]

for j in range(len(sizes)):
    n = sizes[j]
    for i in range(n_of_each):
        for r in R:
            for ins in range(1, 11):
                p = uniform_array(n, 1, r)
                wu = uniform_array(n, 1, r)
                wl = uniform_array(n, 1, r)
                cnt = uniform_array(n, 1, r)
                cl = math.ceil(ins/22*sum(wl)*r/2)
                write_instance(n, p, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl, cnt,
                        {"group":"Count", "class":"uc", "ins":ins/2})
                write_instance(n, wl, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl, cnt,
                        {"group":"Count", "class":"ss-lp", "ins":ins/2})
                write_instance(n, wu, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl, cnt,
                        {"group":"Count", "class":"ss-up", "ins":ins/2})
                write_instance(n, wl, max(max(wl), cl+randint(-10,10)),
                        max(max(wl),cl), wl, wl, cnt,
                        {"group":"Count", "class":"ss-all", "ins":ins/2})
                write_instance(n, p, max(max(wl), cl+randint(-10,10)),
                        max(max(wl),cl), wl, wl, cnt,
                        {"group":"Count", "class":"ss-ul", "ins":ins/2})
