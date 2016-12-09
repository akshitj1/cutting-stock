import numpy as np
from pulp import *
import sys

EPS = 1e-5


def gaussian_elimination(A):
    # only row wise elimination is done
    pivot_idx = -1
    min_ratio = 1e9

    y = A[:, -1]
    x = A[:, -2]
    for idx in range(1, y.shape[0], 1):
        if y[idx] > EPS and x[idx] >= -EPS:
            cur_ratio = x[idx] / y[idx]
            if cur_ratio < min_ratio:
                pivot_idx = idx
                min_ratio = cur_ratio

    if min_ratio < EPS:
        print "min ratio reached 0. add slack vars"
        new_N = x[:] # should be independent of N
        new_N = np.asarray([el if el > EPS else 1 for idx, el in enumerate(new_N)])
        A = np.concatenate((A[:, :-1], new_N.reshape((A.shape[0], 1)), A[:, -1:]), axis=1)
        A = gaussian_elimination(A)
        return A

    if pivot_idx == -1:
        print "Error: pivot was not found"
        sys.exit(1)

    A[pivot_idx, :] /= A[pivot_idx, -1]

    for i in range(A.shape[0]):
        if i != pivot_idx and A[i, -1] > EPS:
            A[i, :] = A[i, :] - A[i, -1] * A[pivot_idx, :]

    # print "matrix after gaussian elimination: \n", A
    return A


def insert_pattern(A):
    B_inv = A[:, 0:A.shape[0]]
    B = np.linalg.inv(B_inv)
    x = A[:, -2]
    y = A[:, -1]
    N_ = B.dot(x)
    P = B.dot(y)

    pivot_idx = -1
    min_ratio = 1e9

    # find the leaving pattern
    for idx in range(1, y.shape[0], 1):
        if y[idx] > EPS and x[idx] >= -EPS:
            cur_ratio = x[idx] / y[idx]
            if cur_ratio < min_ratio:
                pivot_idx = idx
                min_ratio = cur_ratio

    if min_ratio < EPS:
        print "min ratio reached 0. add slack vars"
        # sys.exit(1)

    if pivot_idx == -1:
        print "Error: pivot was not found"
        sys.exit(1)

    B[:, pivot_idx] = P[:]

    # make sure it is invertible

    try:
        B_inv = np.linalg.inv(B)
    except np.linalg.LinAlgError:
        print "Error: matrix is non invertible"
        sys.exit(1)

    N = B_inv.dot(N_)

    # here P is a useless column
    A = np.concatenate((B_inv, N.reshape(-1, 1), P.reshape(-1, 1)), axis=1)
    return A


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


def solve_knapsack(cap, wts, vals, val_lb, geo_idx):
    # todo: replace this ILP approach by dynamic programming

    kp_model = LpProblem("Knapsack problem", LpMaximize)
    lp_vars = ["a-" + str(a) for a in range(len(wts))]
    x = pulp.LpVariable.dicts('assign', lp_vars, lowBound=0, cat=LpInteger)

    # optimization condition
    kp_model += sum([x[lp_var] * vals[idx] for idx, lp_var in enumerate(lp_vars)])

    # knapsack constraint
    for g_idx, geo in enumerate(geo_idx):
        kp_model += sum([ x[lp_vars[idx]] * wts[idx] for idx in range(geo[0], geo[1]+1, 1) ]) <= cap, \
                "do not fill more than capacity for geo %s" % str(g_idx)

    kp_model.writeLP("KnapsackModel.lp")

    kp_model.solve()

    # The status of the solution is printed to the screen
    print "Status:", LpStatus[kp_model.status]

    obj_val = value(kp_model.objective)
    if obj_val - val_lb < EPS:
        return False, None

    print "kp obj val ", obj_val

    pattern = [x[lp_var].value() for lp_var in lp_vars]
    return True, np.asarray(pattern)


def generate_column(stock_lens, stock_costs, opt_coeffs, order_lens, geo_idx):

    for stock_idx in range(len(stock_lens)-1, -1, -1):
        status, pattern = solve_knapsack(stock_lens[stock_idx], order_lens, opt_coeffs, stock_costs[stock_idx], geo_idx)
        if status:
            return pattern, stock_costs[stock_idx]

    # if reached here then it means no better column can be generated
    print("no better column can be generated")
    return None, None


def find_first_greater(el, sorted_list):
    for idx, list_el in enumerate(sorted_list):
        if list_el >= el:
            return idx


def get_duals(patterns, costs, demands):

    cs_model = LpProblem("Cutting stock Master", LpMinimize)
    lp_vars = ["occ-" + str(x) for x in range(patterns.shape[1])]
    x = pulp.LpVariable.dicts('occurrence', lp_vars, lowBound=0, cat=LpContinuous)

    # optimization condition

    cs_model += sum([x[lp_var] * costs[idx] for idx, lp_var in enumerate(lp_vars)])

    # demand constraint
    for order_idx, demand in enumerate(demands):
        cs_model += sum( [x[lp_var] * patterns[order_idx, pattern_idx] for pattern_idx, lp_var in enumerate(lp_vars)] ) == demand, \
        "demand for order %s should be met" % str(order_idx)

    cs_model.writeLP("CuttingStockMaster.lp")

    cs_model.solve()

    obj_val = value(cs_model.objective)
    print ("Z: ", obj_val)

    duals = []
    for name, c in cs_model.constraints.items():
        duals += [c.pi]

    return duals, obj_val


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
    stock_lens = [x for x in range(5, 61, 5)]
    stock_costs = [2*x for x in stock_lens]
    order_lens = [25, 25, 25, 10, 10, 20, 10, 10, 10, 10, 25, 30, 30, 20, 20, 25, 30]
    order_qty = [4, 2, 0, 2, 3, 1, 2, 2, 1, 1,  2, 2, 1, 1, 1, 1, 2]
    geos_idx = [(0, 2), (3, 5), (6, 7), (8, 12), (13, 15), (16, 16)]

    # todo: sort stocks and orders

    geos_idx = geos_idx[:2]
    order_lens = order_lens[:geos_idx[-1][1]+1]
    order_qty = order_qty[:len(order_lens)]

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
    prev_z = 1e9

    for i in range(1000):
        print "-------------------iteration ", i, "--------------------"
        if i == 203:
            break
        B = np.linalg.inv(G[:, 0:G.shape[0]])
        duals, z = get_duals(B[1:, 1:], -1 * B[0, 1:], np.asarray(order_qty))
        if z - prev_z > EPS:
            print "Error: bad column was inserted"
            # break
        prev_z = z
        print "duals are: ", duals

        pattern, pattern_cost = generate_column(stock_lens, stock_costs, np.asarray(duals), order_lens,
                                                geos_idx)

        if pattern is None:
            break

        print "pattern is: ", pattern, " with cost ", pattern_cost

        if pattern_cost > 240:
            print "how????"

        P = np.zeros(N.shape)
        P[0, 0] = -1 * pattern_cost
        P[1:, 0] = pattern
        B_inv = G[:, :G.shape[0]]
        P_ = B_inv.dot(P)
        # print P_

        G[:, -1] = P_[:, 0]
        # print G

        # G = gaussian_elimination(G)
        G = insert_pattern(G)
        # print np.linalg.inv(G[:, 0:G.shape[0]])

    # here we have final matrix G
    # print "final matrix G: \n", G
    B = np.linalg.inv(G[:, 0:G.shape[0]])
    B = np.around(B, decimals=0)
    print "final columns are: \n", B

    get_integral_solution(B[1:, 1:], -1*B[0, 1:], np.asarray(order_qty))

if __name__ == "__main__":
    np.set_printoptions(precision=2)
    gilmore_gomory_column_generation()