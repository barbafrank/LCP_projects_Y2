import numpy as np
import heapq as hq
#import numba

def getDegree(edges):
    if len(edges) == 0:
        return 0
    else:
        degree = 0
        for e in edges:
            degree += e[1]
        return degree

# the function performs a push operation
def push(p, r, c, u, neighbours, neighbours_deg, th):

    # update the vectors (they are passed as pointers)
    delta = r[u]
    p[u] += (1-c)*delta
    r[u] = 0

    # initialize the indicator list that tells
    # approximate pagerank which node needs to be
    # added into the priority queue for the next
    # iteration of the algorithm
    r_above_th = [False]*len(neighbours)

    # update r again and define the values
    # of the indicator list
    for i, e in enumerate(neighbours):
        n = e[0]
        r[n] += c*delta/neighbours_deg[i]
        # a node v must be added if: r(v)/d(v) >= th
        if (neighbours_deg[i] == 0 and r[0] != 0) or r[n]/neighbours_deg[i] >= th:
            r_above_th[i] = True

    return r_above_th

def approximateSimrank(A, N, D, v, c, epsilon, max_iters=200):

    # initialize the approximate PageRank vectors
    p = np.zeros(N)
    r = np.zeros(N)
    r[v] = 1

    # for some reason pyhton decided to create a min pq
    # so we push the cost inverted (in this case d/1 = d)
    pq = []
    hq.heappush(pq, (getDegree(A[v]), v))

    # precompute the inverse of epsilon for efficiency
    th = epsilon/D
    inv_th = D/epsilon

    # iterate over the heap as per defined in [1]
    iters = 0
    while len(pq)>0 and pq[0][0] <= inv_th and iters < max_iters:
        u = hq.heappop(pq)[1]
        neigh = A[u]

        # compute the degrees of the neighbours
        neigh_deg = [0]*len(neigh)
        for i, e in enumerate(neigh):
            n = e[0]
            neigh_deg[i] = getDegree(A[n])

        # call the push function to update p and r and
        # to get the list of neighbours to add to the pq
        r_above_th = push(p, r, c, u, neigh, neigh_deg, th)

        # push the new nodes accordingly to the indicator
        # vector r_above_th returned by push (boolean)
        for i, flag in enumerate(r_above_th):
            if flag:
                n = neigh[i][0]
                # there can never be a math exception, since if r were
                # to be 0 then it would have never been inserted in the list
                hq.heappush(pq, (neigh_deg[i]/r[n], n))

        iters += 1
    # return p
    return p
