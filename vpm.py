# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""
from igraph import *

try:
    from config import *
except ImportError:
    raise ImportError("Config file unavailable")

def sis_vpm_simulate(graph, B, D, c, t):
    """
    Simulates the propagation of a virus with the SIS VPM
    graph: contact network
    B: transmission probability
    D: healing probability
    c: number of initially infected nodes
    t: number of time steps

    The initially infected nodes are chosen from a random uniform
    probability distribution. At each time step, every susceptible (i.e., non-
    infected) node has a 'B' probability of being infected by neighboring
    infected nodes, and every infected node has a 'D' probability of healing
    and becoming susceptible again. The program also calculates the fraction
    of infected nodes at each time step.
    """


if __name__=='__main__':
    g=Graph()
    g.add_vertices(10)
    g.add_edges([(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),
                 (4,5),(5,6),(6,7),(7,8),(8,9),(5,9),(5,8)])
    plot(g)
