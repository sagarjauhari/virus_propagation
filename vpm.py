# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""
from os.path import join
from igraph import Graph
import scipy.linalg

try:
    from config import *
except ImportError:
    raise ImportError("Config file unavailable")

def file2igraph(file):
    """
    Converts graph file into iGraph object, adds artifacts
    """
    with open(file, 'r') as fi:
        v,e = fi.next().split()
        e_list = [(int(i.split()[0]), int(i.split()[1])) for i in list(fi)]
        assert (int(e) == len(e_list)),\
                    "#edges mentioned and # of edges in file differ"
        g = Graph()
        g.add_vertices(int(v))
        g.add_edges(e_list)
        return g

def get_eff_strength(graph, b,d):
    M = len(list(graph.vs))
    eig = scipy.linalg.eigh(graph.get_adjacency().data,
                            eigvals_only=True,
                            eigvals=(M-1, M-1))

    print 'Largest Eigen Value: %0.3f'% (eig)
    return eig*b/d

graph = file2igraph(join(DATA_URL,'static.network'))
eff_strength = get_eff_strength(graph, B1, D1)

