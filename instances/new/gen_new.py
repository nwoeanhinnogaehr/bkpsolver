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
def write_instance(n, p, cu, cl, wu, wl, meta={"group":"New"}):
    global num_inst
    with open("new{}.ki".format(num_inst), "w") as f:
        f.write("{}\n{}\n{}\n".format(n,cl,cu))
        f.write("{}\n".format(" ".join(map(str, wl))))
        f.write("{}\n".format(" ".join(map(str, wu))))
        f.write("{}\n".format(" ".join(map(str, p))))
        for key in meta:
            f.write("{} {}\n".format(key, meta[key]))
        num_inst += 1

sizes = [10, 25, 50, 100, 1000, 10000]
n_of_each = 1
R = [10, 25, 50, 100, 1000]

for j in range(len(sizes)):
    n = sizes[j]
    for i in range(n_of_each):
        for r in R:
            for ins in range(1, 11):
                p = uniform_array(n, 1, r)
                wu = uniform_array(n, 1, r)
                wl = uniform_array(n, 1, r)
                cl = math.ceil(ins/22*sum(wl))
                write_instance(n, p, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl,
                        {"group":"New", "class":"uc", "ins":ins/2})
                write_instance(n, wl, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl,
                        {"group":"New", "class":"ss-lp", "ins":ins/2})
                write_instance(n, wu, max(max(wu), cl+randint(-10,10)),
                        max(max(wl),cl), wu, wl,
                        {"group":"New", "class":"ss-up", "ins":ins/2})
                write_instance(n, wl, max(max(wl), cl+randint(-10,10)),
                        max(max(wl),cl), wl, wl,
                        {"group":"New", "class":"ss-all", "ins":ins/2})
                write_instance(n, p, max(max(wl), cl+randint(-10,10)),
                        max(max(wl),cl), wl, wl,
                        {"group":"New", "class":"ss-ul", "ins":ins/2})
