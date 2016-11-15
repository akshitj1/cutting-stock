from scipy.optimize import linprog
from knapsack import knapsack
import numpy as np


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

if __name__ == "__main__":

    stock = (20, 20)
    # stocks = [Stock(x, x) for x in range(50, 100, 5)]

    orders = [(5, 150), (7, 200), (9, 300)]
    orders = [(length, demand) for length, demand in orders]
    orders = np.array(orders)

    M = len(orders)

    # initial pattern
    patterns = np.zeros((M, M))
    pattern_costs = np.zeros((M))
    for idx, order in enumerate(orders):
        patterns[idx, idx] = stock[0] / order[0]
        pattern_costs[idx] = stock[1]

    solve_master(patterns, orders, stock)


    print patterns

class Stock:
    length = 0
    cost = 0

    def __init__(self, length=0, cost=0):
        self.length = length
        self.cost = cost


class Order:
    length = 0
    demand = 0

    def __init__(self, length=0, demand=0):
        self.length = length
        self.demand = demand



class PatternEntry:
    length = 0
    occurrence = 0

    def __init__(self, length, occurrence):
        self.length = length
        self.occurrence = occurrence


class Pattern:
    # stock used for this pattern
    stock = Stock()
    #items used with their occurrence
    entries = []

    def __init__(self, stock, pattern):
        self.stock = stock
        for item, occ in pattern:
            self.entries.append(PatternEntry(item, occ))


def find_solution(stock, orders):
    # first find initial pattern
    patterns = []
    for idx, order in enumerate(orders):
        pattern = np.array((len(orders)))
        pattern[idx] = PatternEntry(order.length, stock.length/order.length)
