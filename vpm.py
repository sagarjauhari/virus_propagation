# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar
"""
from os.path import join
from igraph import Graph

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

def plot_eig(graph):
    e = np.sort(np.linalg.eig(graph.get_adjacency().data)[0])
    print ["%0.3f" % item for item in e]
    plt.plot(e)
    plt.ylabel('Eigen Values')
    plt.show()

def get_eff_strength(graph, b,d):
    eig = 0 #largest eigen value
    return eig*b/d

a = file2igraph(join(DATA_URL,'static.network'))

