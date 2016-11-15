import numpy as np

def backtrack_DP(values, weights, wt, V, i, j, items):
    # Description - Returns items in the solution using backtracking

    # exit criteria
    if i == 0 or j == 0:
        return items

    # current item not chosen
    if V[i,j] == V[i-1,j]:
        items = backtrack_DP(values, weights, wt, V, i-1, j, items)
    else:
        # current item chosen
        idx = wt.index(wt[j] - weights[i-1])
        items.append(i-1)
        items = backtrack_DP(values, weights, wt, V, i-1, idx, items)

    return items

def knapsack(weights, values, kp):
    modified_weights = []
    modified_values = []
    for idx, weight in enumerate(weights):
        occ = kp / weight
        modified_weights.extend([weight] * occ)
        modified_values.extend([values[idx]] * occ)
    z, pattern = subset_sum(modified_weights, modified_values, kp)
    pattern = [modified_weights[idx] for idx, item in enumerate(pattern) if item > 0]
    zip_pattern = []
    for weight in weights:
        zip_pattern.append(pattern.count(weight))
    return z, zip_pattern


def subset_sum(weights, values, kp):
    # Description - Computes best value for a given capacity using DP

    z = 0
    if kp <= 0:
        return 0, [0] * len(weights)

    # initialize
    pattern = [0] * len(weights)

    wt = range(0,kp+1)
    V = np.zeros((len(weights)+1, len(wt)), dtype = np.double)
    for i in range(1,len(weights)+1):
        for j,w in enumerate(wt):
            if weights[i-1] > w:
                V[i,j] = V[i-1,j]
            else:
                idx = wt.index(w-weights[i-1])
                V[i,j] = max(V[i-1,j], (V[i-1, idx] + values[i-1]))

    r = len(weights)
    c = len(wt)-1
    items = backtrack_DP(values, weights, wt, V, r, c, [])
    for i in items:
        pattern[i] = 1
    z += V[-1,-1]

    return z, pattern
