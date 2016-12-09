"""
cutstock.py:  use gurobi for solving the cutting stock problem.

The instance of the cutting stock problem is represented by the two
lists of m items of size and quantity s=(s_i) and q=(q_i).

The roll size is B.

Given packing patterns t_1, ...,t_k,...t_K where t_k is a vector of
the numbers of items cut from a roll, the problem is reduced to the
following LP:

    minimize   sum_{k} x_k
    subject to sum_{k} t_k(i) x_k >= q_i    for all i
	       x_k >=0			    for all k.

We apply a column generation approch (Gilmore-Gomory approach) in
which we generate cutting patterns by solving a knapsack sub-problem.

Copyright (c) by Joao Pedro PEDROSO and Mikio KUBO, 2010
"""

from gurobipy import *

LOG = True
EPS = 1.e-6


def solveCuttingStock(s, B):
    """solveCuttingStock: use Haessler's heuristic.

    Parameters:
        s - list with item widths
        B - bin capacity

    Returns a solution: list of lists, each of which with the cuts of a roll.
    """
    w = []  # list of different widths (sizes) of items
    q = []  # quantitiy of orders
    for item in sorted(s):
        if w == [] or item != w[-1]:
            w.append(item)
            q.append(1)
        else:
            q[-1] += 1

    t = []  # patterns
    m = len(w)
    # generate initial patterns with one size for each item width
    for i, width in enumerate(w):
        pat = [0] * m  # vector of number of orders to be packed into one roll (bin)
        pat[i] = int(B / width)
        t.append(pat)

    if LOG:
        print "sizes of orders=", w
        print "quantities of orders=", q
        print "roll size=", B
        print "initial patterns", t

    iter = 0
    K = len(t)
    master = Model("LP")  # master LP problem
    x = {}
    for k in range(K):
        x[k] = master.addVar(obj=1, vtype="I", name="x[%d]" % k)
    master.update()

    orders = {}
    for i in range(m):
        coef = [t[k][i] for k in range(K) if t[k][i] > 0]
        var = [x[k] for k in range(K) if t[k][i] > 0]
        orders[i] = master.addConstr(LinExpr(coef, var), ">", q[i], name="Order[%d]" % i)

    master.update()  # must update before calling relax()
    master.Params.OutputFlag = 0  # silent mode
    # master.write("MP" + str(iter) + ".lp")

    while 1:
        iter += 1
        relax = master.relax()
        relax.optimize()
        pi = [c.Pi for c in relax.getConstrs()]  # keep dual variables

        knapsack = Model("KP")  # knapsack sub-problem
        knapsack.ModelSense = -1  # maximize
        y = {}
        for i in range(m):
            y[i] = knapsack.addVar(obj=pi[i], ub=q[i], vtype="I", name="y[%d]" % i)
        knapsack.update()

        L = LinExpr(w, [y[i] for i in range(m)])
        knapsack.addConstr(L, "<", B, name="width")
        knapsack.update()
        # knapsack.write("KP"+str(iter)+".lp")
        knapsack.Params.OutputFlag = 0  # silent mode
        knapsack.optimize()
        if LOG:
            print "objective of knapsack problem:", knapsack.ObjVal
        if knapsack.ObjVal < 1 + EPS:  # break if no more columns
            break

        pat = [int(y[i].X + 0.5) for i in y]  # new pattern
        t.append(pat)
        if LOG:
            print "shadow prices and new pattern:"
            for i, d in enumerate(pi):
                print "\t%5d%12g%7d" % (i, d, pat[i])
            print

        # add new column to the master problem
        col = Column()
        for i in range(m):
            if t[K][i] > 0:
                col.addTerms(t[K][i], orders[i])
        x[K] = master.addVar(obj=1, vtype="I", name="x[%d]" % K, column=col)
        master.update()  # must update before calling relax()
        # master.write("MP" + str(iter) + ".lp")
        K += 1

    # Finally, solve the IP
    if LOG:
        master.Params.OutputFlag = 1  # verbose mode
    master.optimize()

    if LOG:
        print
        print "final solution (integer master problem):  objective =", master.ObjVal
        print "patterns:"
        for k in x:
            if x[k].X > EPS:
                print "pattern", k,
                print "\tsizes:",
                print [w[i] for i in range(m) if t[k][i] > 0 for j in range(t[k][i])],
                print "--> %d rolls" % int(x[k].X + .5)

    rolls = []
    for k in x:
        for j in range(int(x[k].X + .5)):
            rolls.append(sorted([w[i] for i in range(m) if t[k][i] > 0 for j in range(t[k][i])]))
    rolls.sort()
    return rolls


def FFD(s, B):
    """First Fit Decreasing heuristics for the Bin Packing Problem.

    Parameters:
        s - list with item widths
        B - bin capacity
    """
    remain = [B]  # keep list of empty space per bin
    sol = [[]]  # a list ot items (i.e., sizes) on each used bin
    for item in sorted(s, reverse=True):
        for j, free in enumerate(remain):
            if free >= item:
                remain[j] -= item
                sol[j].append(item)
                break
        else:  # does not fit in any bin
            sol.append([item])
            remain.append(B - item)

    return sol


def solveBinPacking(s, B):
    """solveBinPacking: use an IP model to solve the in Packing Problem.

    Parameters:
        s - list with item widths
        B - bin capacity

    Returns a solution: list of lists, each of which with the items in a roll.
    """
    n = len(s)
    U = len(FFD(s, B))  # upper bound of the number of bins
    model = Model("bpp")
    # setParam("MIPFocus",1)
    x = {}
    listX = {}  # for SOS constraints	!!!!!!!!!!!!!!!!!!!
    for i in range(n):
        listX[i] = []
        for j in range(U):
            x[i, j] = model.addVar(vtype="B", name="x[%d,%d]" % (i, j))
            listX[i].append(x[i, j])
    y = {}
    for j in range(U):
        y[j] = model.addVar(obj=1, vtype="B", name="y[%d]" % j)
    model.update()

    # assignment constraints
    for i in range(n):
        var = [x[i, j] for j in range(U)]
        coef = [1] * U
        model.addConstr(LinExpr(coef, var), "=", 1, name="cnstr1[%d]" % i)

    # bin capacity constraints
    for j in range(U):
        var = [x[i, j] for i in range(n)]
        coef = [s[i] for i in range(n)]
        model.addConstr(LinExpr(coef, var), "<", LinExpr(B, y[j]), name="cnstr2[%d]" % j)

    # tighten assignment constraints
    for j in range(U):
        for i in range(n):
            model.addConstr(x[i, j], "<", y[j], name="cnstr3[%d,%d]" % (i, j))

            ##    # tie breaking constraints
            ##    for j in range(U-1):
            ##        lin=LinExpr()
            ##        lin.addTerms(1,y[j])
            ##        lin.addTerms(-1,y[j+1])
            ##        model.addConstr(lhs=lin,sense=GRB.GREATER_EQUAL,rhs=0,name="constraint4_"+str(j))
            ##
            ##    # SOS constraints
            ##    for i in range(n):
            ##        model.addSOS(1,listX[i])

    if not LOG:
        model.Params.outputflag = 0
    model.optimize()

    bins = [[] for i in range(U)]
    for (i, j) in x:
        if x[i, j].X > EPS:
            bins[j].append(s[i])
    for i in range(bins.count([])):
        bins.remove([])
    for b in bins:
        b.sort()
    bins.sort()
    return bins


def CuttingStockExample1():
    """CuttingStockExample1: create toy instance for the cutting stock problem."""
    B = 110  # roll width (bin size)
    w = [20, 45, 50, 55, 75]  # width (size) of orders (items)
    q = [48, 35, 24, 10, 8]  # quantitiy of orders
    s = []
    for j in range(len(w)):
        for i in range(q[j]):
            s.append(w[j])
    return s, B


def CuttingStockExample2():
    """CuttingStockExample2: create toy instance for the cutting stock problem."""
    B = 9  # roll width (bin size)
    w = [2, 3, 4, 5, 6, 7, 8]  # width (size) of orders (items)
    q = [4, 2, 6, 6, 2, 2, 2]  # quantitiy of orders
    s = []
    for j in range(len(w)):
        for i in range(q[j]):
            s.append(w[j])
    return s, B


def DiscreteUniform(n=10, LB=1, UB=99, B=100):
    """DiscreteUniform: create random, uniform instance for the bin packing problem."""
    import random
    random.seed(1)
    B = 100
    s = [0] * n
    for i in range(n):
        s[i] = random.randint(LB, UB)
    return s, B


if __name__ == "__main__":
    # s,B = CuttingStockExample1()
    s, B = CuttingStockExample2()
    # n = 500
    # B = 100
    # s,B = DiscreteUniform(n,18,50,B)

    ffd = FFD(s, B)
    print "\n\n\nSolution of FFD:"
    print ffd
    print len(ffd), "bins"

    print "\n\n\nCutting stock problem:"
    rolls = solveCuttingStock(s, B)
    print len(rolls), "rolls:"
    print rolls

    print "\n\n\nBin packing problem:"
    bins = solveBinPacking(s, B)
    print len(bins), "bins:"
    print bins