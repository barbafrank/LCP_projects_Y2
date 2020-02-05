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
    cdef dict idxs_map
    cdef list raw_edge, A, original_idxs
    cdef edge_t edge
    cdef (int, double) clean_edge
    cdef int N, n

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
        idxs_map = {n:i for i, n in enumerate(edge_dict.keys())}
        original_idxs = list(idxs_map.keys())
        N = len(edge_dict.keys())
        A = [[]]*N

        for n in range(N):
            A[n] = []
            if original_idxs[n] in edge_dict:
                A[n] = [(idxs_map[clean_edge[0]], clean_edge[1]) for clean_edge in edge_dict[original_idxs[n]]]

        return (A, N, original_idxs)
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
        idxs_map = {n:i for i, n in enumerate(edge_dict.keys())}
        original_idxs = list(idxs_map.keys())
        N = len(edge_dict.keys())
        A = [[]]*N

        for n in range(N):
            A[n] = []
            if original_idxs[n] in edge_dict:
                A[n] = [(idxs_map[clean_edge[0]], clean_edge[1]) for clean_edge in edge_dict[original_idxs[n]]]

        return (A, N, original_idxs)
    else:
        raise Exception("No such encoding type:", enc_type)
