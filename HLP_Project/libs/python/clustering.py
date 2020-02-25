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

def cluster1(L_norm, centroids, alpha):
    # initialize the output vector with the clusters
    clusters = - np.ones(len(L_norm))
    
    # initialize the heap queues related to the clusters
    queues = []
    
    # insert centroids in the queues, e.g. (expansions, (node, rank))
    for i in range(centroids.size):
        queue = []
        hq.heappush(queue, (0, (centroids[i], 1)))
        queues.append(queue)
        
        # set the centroid cluster
        clusters[centroids[i]] = centroids[i]
    
    # iterate until all the nodes are assigned to the clusters   
    exit_ext = False
    while not exit_ext:
        cumulative_length = 0
        for queue in queues:
            # assign the new expansions
            for i, node in enumerate(queue):
                # ensure at least one expansion per node
                expansions = max(1, np.round(alpha*node[1][1])) # TODO: eventually change the round function
                
                # modify the queue element with the new expansions
                queue[i] = (
                    - expansions,
                    (node[1][0], node[1][1])
                )
            
            # restore the queue
            hq.heapify(queue)
        
            # expand the cluster boundary toward the high ranked nodes
            exit_int = False
            while (not exit_int) and len(queue)>0:
                # pop a node that should be expanded 
                node = hq.heappop(queue)
                
                # check the termination condition
                expansions_opp = node[0]
                if expansions_opp==0:
                    # no more nodes to expand
                    exit_int = True
                else:
                    # reduce the remained number of expansions
                    expansions_opp += 1
                    
                    # expand the neighbourhood and assign the new expansions
                    for edge in L_norm[node[1][0]]:
                        #print(neighbour)
                        neighbour = edge[0]
                        # assign the node to the cluster
                        if clusters[neighbour]<0:
                            #print("assigned")
                            clusters[neighbour] = clusters[node[1][0]]
                            
                            # add the node to the queue
                            rank = edge[1]
                            hq.heappush(
                                queue,
                                (
                                    np.round(expansions_opp*rank),
                                    (neighbour, rank)
                                )
                            )
            cumulative_length += len(queue) 
        exit_ext = cumulative_length==0
    
    # return the vector with the clusters
    return clusters

def cluster3(L):
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