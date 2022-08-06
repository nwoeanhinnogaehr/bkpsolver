import random
n = 15
for i in range(100):
    with open("genssa15_{}.ki".format(i), "w") as f:
        f.write(str(n) + "\n")
        f.write(str(random.randint(100, 500)) + "\n")
        f.write(str(random.randint(100, 500)) + "\n")
        s = [random.randint(1, 100) for _ in range(n)]
        for j in range(n):
            f.write(str(s[j]) + " ")
        f.write("\n")
        for j in range(n):
            f.write(str(s[j]) + " ")
        f.write("\n")
        for j in range(n):
            f.write(str(s[j]) + " ")
        f.write("\n")
