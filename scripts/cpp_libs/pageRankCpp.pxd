#!python
#distutils: language=c++
#cython: language_level=3

from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool as cbool

cdef extern from "page_rank.cpp":
    pass

cdef extern from "page_rank.h":
    ctypedef pair[int, double] edge_t
    ctypedef vector[edge_t] nodelist_t
    ctypedef vector[nodelist_t] edgelist_t

    double getDegree(nodelist_t edges)
    edgelist_t localPageRank(edgelist_t A, double c, double epsilon, int max_iter, cbool return_only_neighbours)
