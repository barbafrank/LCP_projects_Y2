import numpy as np

# Method to convert the edgelist back to a dense matrix (Just for debugging purposes)
def list2matrix(edgelist):
    N = len(edgelist)
    A_mat = np.zeros((N,N))
    for i, node in enumerate(edgelist):
        for neighbour in node:
            A_mat[i, neighbour[0]] = neighbour[1]
    return A_mat

def getInOutDegree(edgelist):
    N = len(edgelist)
    outDegs = [0.]*N
    inDegs = [0.]*N

    # for each node and list
    for i, n in enumerate(edgelist):
        # initialize out degree
        out_deg = 0.
        # for each weight
        for edge in n:
            # update out degree
            out_deg += edge[1]
            # and in degree
            inDegs[edge[0]] += edge[1]
        # set out degree to the correct value
        outDegs[i] = out_deg

    return inDegs, outDegs

# compute the transposition of the edgelist
def transposeEdgelist(edgelist):
    N = len(edgelist)
    edgelist_t = [[]]*N

    for i in range(len(edgelist_t)):
        edgelist_t[i] = []

    for i, nodelist in enumerate(edgelist):
        for edge in nodelist:
            edgelist_t[edge[0]].append((i, edge[1]))

    return edgelist_t

def reorderEdgelist(edgelist, new_idxs):
    N = len(edgelist)
    edgelist_t = [[]]*N

    for i in range(len(edgelist_t)):
        edgelist_t[i] = []

    for i, nodelist in enumerate(edgelist):
        for edge in nodelist:
            edgelist_t[new_idxs[i]].append((new_idxs[edge[0]], edge[1]))

    return edgelist_t
