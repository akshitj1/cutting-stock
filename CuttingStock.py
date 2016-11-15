import numpy as np
from scipy.optimize import linprog
from knapsack import knapsack
from pulp import *

# find entering column
def find_entering_pattern(patterns, pattern_costs, stock, orders):
    A = patterns
    C = pattern_costs

    # knapsack maximization
    c = stock[1]
    B = np.matrix(C) * np.linalg.inv(np.matrix(A))

    # knapsack constraint
    l = stock[0]
    L = [order[0] for order in orders]

    z, pattern = knapsack(L, B.tolist()[0], l)

    # our new pattern should give greater c where c is the stock cost
    print "new c is: ", z
    print "find leaving pattern from \n", A, " for entering pattern ", pattern
    return pattern

def insert_new_pattern(patterns, entering_pattern, pattern_occ):
    min_ratio = float('inf')
    leaving_pattern_idx = -1
    num_matrix = np.linalg.inv(np.matrix(patterns)) * np.matrix(entering_pattern).T
    print num_matrix
    for idx, occ in enumerate(entering_pattern):
        if occ <= 0:
            continue
        cur_ratio = num_matrix[idx, 0] / occ
        if min_ratio > cur_ratio:
            min_ratio = cur_ratio
            leaving_pattern_idx = idx

    if leaving_pattern_idx == -1:
        print "no leaving pattern found"

    patterns[:, leaving_pattern_idx] = entering_pattern[:]
    return patterns


def solve_master(patterns, orders, stock):
    M = len(orders)

    opt_coeff = np.zeros((M))
    opt_coeff[:] = stock[1]

    # patterns
    A_ub = np.zeros((M, M))
    # demand
    b_ub = np.zeros((M))

    for idx in range(M):
        b_ub[idx] = orders[idx, 1]

    b_ub = -1 * b_ub

    for loop_count in range(0, 2):
        print "executing loop ", loop_count, "----------------------------------------------------"
        A_ub = -1*patterns

        res = linprog(opt_coeff, A_ub=A_ub, b_ub=b_ub)

        print res
        print "check if demand is fulfilled ", np.matrix(patterns) * np.matrix(res.x).T
        patterns_occ = res.x
        min_cost = res.fun

        entering_pattern = find_entering_pattern(patterns, pattern_costs, stock, orders)
        patterns = insert_new_pattern(patterns, entering_pattern, patterns_occ)
        print "new pattern matrix is \n", patterns
#
# if __name__ == "__main__":
#     gilmore_gomory_column_generation()
#
#     stock = (20, 20)
#     # stocks = [Stock(x, x) for x in range(50, 100, 5)]
#
#     orders = [(5, 150), (7, 200), (9, 300)]
#     orders = [(length, demand) for length, demand in orders]
#     orders = np.array(orders)
#
#     M = len(orders)
#
#     # initial pattern
#     patterns = np.zeros((M, M))
#     pattern_costs = np.zeros((M))
#     for idx, order in enumerate(orders):
#         patterns[idx, idx] = stock[0] / order[0]
#         pattern_costs[idx] = stock[1]
#
#     solve_master(patterns, orders, stock)


    # print patterns

# class Stock:
#     length = 0
#     cost = 0
#
#     def __init__(self, length=0, cost=0):
#         self.length = length
#         self.cost = cost
#
#
# class Order:
#     length = 0
#     demand = 0
#
#     def __init__(self, length=0, demand=0):
#         self.length = length
#         self.demand = demand
#
#
#
# class PatternEntry:
#     length = 0
#     occurrence = 0
#
#     def __init__(self, length, occurrence):
#         self.length = length
#         self.occurrence = occurrence
#
#
# class Pattern:
#     # stock used for this pattern
#     stock = Stock()
#     #items used with their occurrence
#     entries = []
#
#     def __init__(self, stock, pattern):
#         self.stock = stock
#         for item, occ in pattern:
#             self.entries.append(PatternEntry(item, occ))
#

# def find_solution(stock, orders):
#     # first find initial pattern
#     patterns = []
#     for idx, order in enumerate(orders):
#         pattern = np.array((len(orders)))
#         pattern[idx] = PatternEntry(order.length, stock.length/order.length)

EPS = 1e-5

def gaussian_elimination(A):
    # only row wise elimination is done
    pivot_idx = -1
    for i in range(A.shape[0]):
        if i!=0 and A[i, -1] > EPS:
            if pivot_idx == -1:
                pivot_idx = i
            else:
                if A[pivot_idx, -1] > A[i, -1]:
                    pivot_idx = i

    A[pivot_idx, :] /= A[pivot_idx, -1]
    # print A

    for i in range(A.shape[0]):
        if i != pivot_idx and A[i, -1] > EPS:
            A[i, :] = A[i, :] - A[i, -1] * A[pivot_idx, :]

    print "matrix after gaussian elimination: \n", A


def solve_adhoc(cap, wts, vals, val_lb):
    # first try adhoc
    # sort in descending ratio
    ratios = [(vals[idx]/wts[idx], idx) for idx, val in enumerate(vals)]
    ratios = sorted(ratios, key=lambda ratio: ratio[0], reverse=True)

    pattern = np.zeros(len(vals))

    rem_capacity=cap
    cur_val = 0
    for idx, v in enumerate(ratios):
        if rem_capacity <=0:
            break
        item_idx = v[1]
        occ = rem_capacity / wts[item_idx]
        pattern[item_idx] = occ
        cur_val += occ * vals[item_idx]
        rem_capacity -= occ * wts[item_idx]

    if cur_val > val_lb:
        return True, pattern
    else:
        return False, pattern


def solve_knapsack(cap, wts, vals, val_lb):
    # todo: replace this ILP approach by dynamic programming

    kp_model = LpProblem("Knapsack problem", LpMaximize)
    lp_vars = ["a-" + str(a) for a in range(len(wts))]
    x = pulp.LpVariable.dicts('assign', lp_vars, lowBound=0, cat=LpInteger)

    # optimization condition
    kp_model += sum([x[lp_var] * vals[idx] for idx, lp_var in enumerate(lp_vars)])

    # knapsack constraint
    kp_model += sum([x[lp_var] * wts[idx] for idx, lp_var in enumerate(lp_vars)]) <= cap, \
                "do not fill more than capacity"

    kp_model.writeLP("KnapsackModel.lp")

    kp_model.solve()

    # The status of the solution is printed to the screen
    print "Status:", LpStatus[kp_model.status]

    obj_val = value(kp_model.objective)
    if obj_val <= val_lb:
        return False, None

    pattern = [ x[lp_var].value() for lp_var in lp_vars ]
    return True, np.asarray(pattern)


def generate_column(stock_lens, stock_costs, opt_coeffs, order_lens):
    for stock_idx in range(len(stock_lens)-1, -1, -1):
        status, pattern =solve_adhoc(stock_lens[stock_idx], order_lens, opt_coeffs, stock_costs[stock_idx])
        if status:
            return pattern, stock_costs[stock_idx]

    # if reached here then it means ad hoc gave no solution
    print("no solution was found using ad hoc method. trying knapsack now.")

    for stock_idx in range(len(stock_lens)-1, -1, -1):
        status, pattern =solve_knapsack(stock_lens[stock_idx], order_lens, opt_coeffs, stock_costs[stock_idx])
        if status:
            return pattern, stock_costs[stock_idx]

    # if reached here then it means no better column can be generated
    print("no better column can be generated")
    return None, None


def find_first_greater(el, sorted_list):
    for idx, list_el in enumerate(sorted_list):
        if list_el>=el:
            return idx


def get_integral_solution(patterns, costs, demands):

    cs_model = LpProblem("Cutting stock problem", LpMinimize)
    lp_vars = ["occ-" + str(x) for x in range(patterns.shape[1])]
    x = pulp.LpVariable.dicts('occurrence', lp_vars, lowBound=0, cat=LpInteger)

    # optimization condition

    cs_model += sum([x[lp_var] * costs[idx] for idx, lp_var in enumerate(lp_vars)])

    # demand constraint
    for order_idx, demand in enumerate(demands):
        cs_model += sum( [x[lp_var] * patterns[order_idx, pattern_idx] for pattern_idx, lp_var in enumerate(lp_vars)] ) >= demand, \
        "demand for order %s should be met" % str(order_idx)

    cs_model.writeLP("CuttingStockModel.lp")

    cs_model.solve()

    # The status of the solution is printed to the screen
    print "Status:", LpStatus[cs_model.status]

    obj_val = value(cs_model.objective)
    print ("demand was met with cost ", obj_val)

    pattern_occ = [x[lp_var].value() for lp_var in lp_vars]
    print "pattern occurences ", pattern_occ


def gilmore_gomory_column_generation():
    stock_lens = [5, 6, 9]
    stock_costs = [x+1 for x in stock_lens]
    order_lens = [2, 3, 4]
    order_qty = [20, 10, 20]

    # todo: sort stocks and orders

    M = len(order_lens)

    B = np.zeros((M+1, M+1))

    B[0, 0] = 1

    N_ = np.zeros((M + 1, 1))

    # fill first row with stock used costs
    for order_idx in range(M):
        order_len = order_lens[order_idx]
        # find just greater stock
        stock_idx = find_first_greater(order_len, stock_lens)
        stock_len = stock_lens[stock_idx]
        stock_cost = stock_costs[stock_idx]
        B[0, order_idx+1] = -1 * stock_cost
        B[order_idx+1, order_idx+1] = stock_len / order_len

        N_[order_idx+1] = order_qty[order_idx]

    print (B)
    print (N_)

    B_inv = np.linalg.inv(B)
    print(B_inv)

    N = B_inv.dot(N_)
    print N

    G = np.concatenate((B_inv, N, N), axis=1)

    for i in range(100):
        print "-------------------iteration ", i, "--------------------"

        pattern, pattern_cost = generate_column(stock_lens, stock_costs, G[0, 1 : 1 + len(order_lens)], order_lens)

        if pattern is None:
            break

        print "pattern is: ", pattern, " with cost ", pattern_cost

        P = np.zeros(N.shape)
        P[0, 0] = -1 * pattern_cost
        P[1:, 0] = pattern
        B_inv = G[:, :G.shape[0]]
        P_ = B_inv.dot(P)
        print P_

        G[:, -1] = P_[:, 0]
        print G

        gaussian_elimination(G)

    # here we have final matrix G
    print "final matrix G: \n", G
    B = np.linalg.inv(G[:, 0:G.shape[0]])
    B = np.around(B, decimals=0)
    print "final columns are: \n", B

    get_integral_solution(B[1:, 1:], -1*B[0, 1:], np.asarray(order_qty))

if __name__ == "__main__":
    np.set_printoptions(precision=2)
    gilmore_gomory_column_generation()