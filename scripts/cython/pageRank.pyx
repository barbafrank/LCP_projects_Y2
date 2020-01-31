#!python
#distutils: language=c++
#cython: language_level=3, boundscheck=False

import numpy as np
cimport numpy as np
import heapq as hq
from scripts.cython.utils import getDegree

#******************************************************************************
# other pagerank stuffs
# the function performs a push operation
def push_other(np.ndarray p, np.ndarray r, double c, int u, list neighbours, np.ndarray neighbours_deg, double th):

    cdef int n
    cdef int i
    cdef (int, double) e

    # update the vectors (they are passed as pointers)
    cdef double delta = r[u]
    p[u] += (1-c)*delta
    r[u] = 0

    # initialize the indicator list that tells
    # approximate pagerank which node needs to be
    # added into the priority queue for the next
    # iteration of the algorithm
    cdef np.ndarray r_above_th = np.zeros([len(neighbours)], dtype=np.bool_)

    # update r again and define the values
    # of the indicator list
    for i, e in enumerate(neighbours):
        n = e[0]
        r[n] += c*delta/neighbours_deg[i]
        # a node v must be added if: r(v)/d(v) >= th
        if (neighbours_deg[i] == 0 and r[0] != 0) or r[n]/neighbours_deg[i] >= th:
            r_above_th[i] = True

    return list(r_above_th)

def approximateSimrank_other(list A, int N, double D, int v, double c, double epsilon, int max_iters=200):

    # initialize the approximate PageRank vectors
    cdef np.ndarray p = np.zeros([N], dtype=np.float64)
    cdef np.ndarray r = np.zeros([N], dtype=np.float64)

    cdef list pq = []
    cdef double th
    cdef double inv_th
    cdef int iters = 0
    cdef int u, i, n
    cdef list neigh
    cdef np.ndarray neigh_deg
    cdef (int, double) e
    cdef list r_above_th
    cdef int flag

    r[v] = 1

    # for some reason pyhton decided to create a min pq
    # so we push the cost inverted (in this case d/1 = d)
    hq.heappush(pq, (getDegree(A[v]), v))

    # precompute the inverse of epsilon for efficiency
    th = epsilon/D
    inv_th = D/epsilon

    # iterate over the heap as per defined in [1]
    while len(pq)>0 and pq[0][0] <= inv_th and iters < max_iters:
        u = hq.heappop(pq)[1]
        neigh = A[u]

        # compute the degrees of the neighbours
        neigh_deg = np.zeros([len(neigh)], dtype=np.float64)
        for i, e in enumerate(neigh):
            n = e[0]
            neigh_deg[i] = getDegree(A[n])

        # call the push function to update p and r and
        # to get the list of neighbours to add to the pq
        r_above_th = push_other(p, r, c, u, neigh, neigh_deg, th)

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

# define the lazy pagerank
#*********************************************************+********************
# the function performs a push operation
def push(np.ndarray p, np.ndarray r, double alpha, int u, double du,
                    list neighbours, np.ndarray neighbours_deg, double epsilon):

    cdef int n
    cdef int i
    cdef (int, double) e

    # update the vectors (they are passed as pointers)
    p[u] += alpha*r[u]
    r[u] *= (1.0 - alpha)/2.

    # initialize the indicator list that tells
    # approximate pagerank which node needs to be
    # added into the priority queue for the next
    # iteration of the algorithm
    cdef np.ndarray r_above_th = np.zeros([len(neighbours)], dtype=np.bool_)

    # update r again and define the values
    # of the indicator list
    for i, e in enumerate(neighbours):
        n = e[0]
        r[n] += (1.0 - alpha)*r[u]/(2.0*du)
        # a node v must be added if: r(v)/d(v) >= epsilon
        if (neighbours_deg[i] == 0 and r[n] != 0) or r[n]/neighbours_deg[i] >= epsilon:
            r_above_th[i] = True

    return list(r_above_th)

def approximateSimrank(list A, int v, double alpha, double epsilon, int max_iters=200,
                            int return_residual=False, int return_only_neighbours=False):
    # compute N
    cdef int N = len(A)

    # initialize the approximate PageRank vectors
    cdef np.ndarray p = np.zeros([N], dtype=np.float64)
    cdef np.ndarray r = np.zeros([N], dtype=np.float64)

    cdef list pq = []
    cdef double inv_epsilon, du
    cdef list v_neighs, neigh
    cdef int iters
    cdef int u, i, flag
    cdef np.ndarray neigh_deg
    cdef (int, double) e
    cdef list r_above_th

    r[v] = 1

    # for some reason pyhton decided to create a min pq
    # so we push the cost inverted (in this case d/1 = d)
    hq.heappush(pq, (getDegree(A[v]), v))

    # precompute the inverse of epsilon for efficiency
    inv_epsilon = 1/epsilon

    # extract the original neighbours
    v_neighs = [n[0] for n in A[v]]

    # iterate over the heap as per defined in [1]
    iters = 0
    while len(pq)>0 and pq[0][0] <= inv_epsilon and iters<max_iters:
        u = hq.heappop(pq)[1]
        neigh = A[u]

        # compute the degrees of the neighbours
        neigh_deg = np.zeros([len(neigh)], dtype=np.float64)
        du = 0
        for i, e in enumerate(neigh):
            n = e[0]
            du += e[1]
            neigh_deg[i] = getDegree(A[n])

        # call the push function to update p and r and
        # to get the list of neighbours to add to the pq
        r_above_th = push(p, r, alpha, u, du, neigh, neigh_deg, epsilon)

        # push the new nodes accordingly to the indicator
        # vector r_above_th returned by push (boolean)
        for i, flag in enumerate(r_above_th):
            if flag:
                n = neigh[i][0]
                # there can never be a math exception, since if r were
                # to be 0 then it would have never been inserted in the list
                hq.heappush(pq, (neigh_deg[i]/r[n], n))

        iters += 1

    if return_only_neighbours:
        p = p[v_neighs]
        r = r[v_neighs]

    if not return_residual:
        return p
    else:
        return p,r


# the function creates the L matrix iterating locally the approximate simrank
#def localPageRank(A, N, D, c, epsilon=1e-5, max_iters=200):
def localPageRank(list A, double c, double epsilon=1e-5, int max_iters=200,
                                                int return_only_neighbours=False):

    cdef int N = len(A)
    #L = np.zeros((N, N))
    cdef list L = [None]*N

    # equivalent lazy pagerank constant
    cdef double alpha = 2*c/(1+c)
    # andersen's paper inverts alpha
    alpha = 1-alpha

    cdef int i, k
    cdef list node
    cdef np.ndarray p
    cdef (int, double) neighbour

    for i, node in enumerate(A):
        #p = approximateSimrank(A, N, D, i, c, epsilon, max_iters)
        p = approximateSimrank(A, i, alpha, epsilon, max_iters, False,
                                                            return_only_neighbours)

        L[i] = []

        for k, neighbour in enumerate(node):
            #L[i,neighbour[0]] = p[neighbour[0]]
            if return_only_neighbours:
                L[i].append((neighbour[0], p[k]))
            else:
                L[i].append((neighbour[0], p[neighbour[0]]))

    return L
