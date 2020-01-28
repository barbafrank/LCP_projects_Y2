import numpy as np

# define the routine that parses a text file
# into the efficient adjacency representation
# needed by the other functions
def edgelistParser(path, enc_type="list"):

    # initialize the edge dictionary
    edge_dict = {}

    # exprects the dump of the upprer tringular part of
    # a symmetric adjacency matrix, meaning the network
    # must be undirected
    if enc_type == "list":
        # the edgelist is a text file where each row
        # is a link between nodes i and j
        # if there is no third number then the matrix is
        # assumed unweighted

        # open the file
        with open(path, "r") as f:
            # loop over its lines
            for line in f:
                # split into chunks
                raw_edge = line.split()
                assert len(raw_edge)>1, "There must be at least 2 numbers per line"
                # check if weighted
                if len(raw_edge)>2:
                    edge = tuple([int(val) for val in raw_edge[:-1]])
                else:
                    edge = tuple([int(val) for val in raw_edge]) + (1,)

                # add the edges to the dictionary
                if edge[0] in edge_dict:
                    edge_dict[edge[0]].append((edge[1], edge[2]))
                else:
                    edge_dict[edge[0]] = [(edge[1], edge[2])]
                if edge[0] != edge[1]:
                    if edge[1] in edge_dict:
                        edge_dict[edge[1]].append((edge[0], edge[2]))
                    else:
                        edge_dict[edge[1]] = [(edge[0], edge[2])]

        # finally convert to list and return
        N = max(edge_dict.keys())+1
        A = [[]]*N
        D = 0

        for n in range(N):
            if n in edge_dict:
                A[n] = edge_dict[n]
                D += len(A[n])

        return (A, N, D)
    # expects duplicated links, as per a undirected
    # network saved as directed
    elif enc_type == "raw_list":
        # the edgelist is a text file where each row
        # is a link between nodes i and j
        # if there is no third number then the matrix is
        # assumed unweighted

        # open the file
        with open(path, "r") as f:
            # loop over its lines
            for line in f:
                # split into chunks
                raw_edge = line.split()
                assert len(raw_edge)>1, "There must be at least 2 numbers per line"
                # check if weighted
                if len(raw_edge)>2:
                    edge = tuple([int(val) for val in raw_edge[:-1]])
                else:
                    edge = tuple([int(val) for val in raw_edge]) + (1,)

                # add the edges to the dictionary
                if edge[0] in edge_dict:
                    edge_dict[edge[0]].append((edge[1], edge[2]))
                else:
                    edge_dict[edge[0]] = [(edge[1], edge[2])]

        # finally convert to list and return
        N = max(edge_dict.keys())+1
        A = [[]]*N
        D = 0

        for n in range(N):
            if n in edge_dict:
                A[n] = edge_dict[n]
                D += len(A[n])

        return (A, N, D)
    else:
        raise Exception("No such encoding type:", enc_type)

# Method to convert the edgelist back to a dense matrix (Just for debugging purposes)
def list2matrix(edgelist):
    count = 0
    N = len(edgelist)
    A_mat = np.zeros((N,N))
    for i,node in enumerate(edgelist):
        for neighbour in node:
            if A_mat[i, neighbour[0]] != 0:
                count += 1;
            A_mat[i, neighbour[0]] = neighbour[1]
    print(count)
    return A_mat
