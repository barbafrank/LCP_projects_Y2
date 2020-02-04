#!python
#distutils: language=c++
#cython: language_level=3, boundscheck=False
from scripts.cython.utils import getDegree

ctypedef (int, int, double) edge_t

# define the routine that parses a text file
# into the efficient adjacency representation
# needed by the other functions
def edgelistParser(str path, str enc_type="list"):

    # initialize the edge dictionary
    cdef dict edge_dict = {}
    cdef list raw_edge, A
    cdef edge_t edge
    cdef (int, double) clean_edge
    cdef int N, D, min_n

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
                    edge = tuple([int(val) for val in raw_edge[:-2]]) + (float(raw_edge[-2]),)
                else:
                    edge = tuple([int(val) for val in raw_edge]) + (1.,)

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
        min_n = min(edge_dict.keys())
        N = max(edge_dict.keys())-min_n+1
        A = [[]]*N
        D = 0

        for n in range(N):
            if n+min_n in edge_dict:
                A[n] = [(clean_edge[0]-min_n, clean_edge[1]) for clean_edge in edge_dict[n+min_n]]
                D += getDegree(A[n])

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
        min_n = min(edge_dict.keys())
        N = max(edge_dict.keys())-min_n+1
        A = [[]]*N
        D = 0

        for n in range(N):
            if n+min_n in edge_dict:
                A[n] = [(clean_edge[0]-min_n, clean_edge[1]) for clean_edge in edge_dict[n+min_n]]
                D += getDegree(A[n])

        return (A, N, D)
    else:
        raise Exception("No such encoding type:", enc_type)
