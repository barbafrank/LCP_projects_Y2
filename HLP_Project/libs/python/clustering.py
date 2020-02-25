import numpy as np
import pandas as pd
import networkx as nx
import heapq as hq

def generate_edges(positions, radius=1, weighted=False):
    # initialize the edgelist
    N = positions.shape[0]
    idxs = np.arange(N)
    A = [[]]*N
    
    # generate the edgelist
    for i in range(N):
        point = positions[i,:]
        dists = np.linalg.norm(point-positions, axis=1)
        neighs = idxs[np.logical_and(dists<=radius, idxs!=i)]
        A[i] = []
        for j, neigh in enumerate(neighs):
            if weighted:
                A[i].append((neigh, 1 + 1./dists[idxs!=i][j]))
            else:   
                A[i].append(neigh)

    # return the generated edgelist
    return A

def build_network(path, sep='\t', radius=1, weighted=False):
    df = pd.read_csv(path, sep=sep, names=["x", "y", "cluster"])
    positions = np.array(df[["x", "y"]])
    A = generate_edges(positions, radius, weighted)
    return A, positions

def cluster4(L):
    clusters = -np.ones(len(L), dtype=int)
    parents = np.zeros(len(L), dtype=int)
    
    # for each node and its neighbours
    for idx in range(len(L)):
        # if it was not already assigned to a cluster:
        if clusters[idx] == -1:
            clusters[idx] = idx
            exit_flag = False
            curr_node = idx
            if len(L[curr_node]) == 0:
                exit_flag = True
            while not exit_flag:
                next_node = max(L[curr_node], key=lambda x:x[1])[0]
                if clusters[next_node] == -1:
                    clusters[next_node] = idx
                    parents[next_node] = curr_node
                    curr_node = next_node
                else:
                    # if the path didn't close on itself
                    if not clusters[next_node] == idx:
                        # walk backwards and reassign the nodes to the correct cluster
                        while curr_node != idx:
                            clusters[curr_node] = clusters[next_node]
                            curr_node = parents[curr_node]
                        clusters[idx] = clusters[next_node]
                    exit_flag = True

    # return the clusters
    return clusters