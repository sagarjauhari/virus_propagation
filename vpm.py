# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""
from igraph import *
import random

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
    c: fraction of initially infected nodes
    t: number of time steps

    The initially infected nodes are chosen from a random uniform
    probability distribution. At each time step, every susceptible (i.e., non-
    infected) node has a 'B' probability of being infected by neighboring
    infected nodes, and every infected node has a 'D' probability of healing
    and becoming susceptible again. The program also calculates the fraction
    of infected nodes at each time step.

    During each time interval t, an infected node i tries to infect its
    neighbors with probability B. At the same time, i may be cured with
    probability D.
    """
    assert 0<=c<=1, ' c should be between 0 and 1'
    N = len(list(graph.vs)) # Number of nodes
    infected = set(random.sample(xrange(N),int(c*N)))

    #Start simulation
    for _i in range(t):
        # For each infected node, infect its neighbors with probability B
        infected_new = set()
        for n in infected:
            nbrs = graph.neighbors(n)
            nbrs_infected = random.sample(nbrs, int(B*len(nbrs)))
            for j in nbrs_infected:
                infected_new.add(j)

        # For all the nodes which were in infected state before this time step:
        # cure them with probability D
        cured = random.sample(infected, int(math.ceil(D*len(infected))))

        # So, now we have infected (old), infected_new & cured
        for n in infected_new:
            infected.add(n)

        for n in cured:
            infected.remove(n)

        print 'Infected Nodes: %d'%(len(infected))




if __name__=='__main__':
    g=Graph()
    g.add_vertices(10)
    g.add_edges([(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),
                 (4,5),(5,6),(6,7),(7,8),(8,9),(5,9),(5,8)])
    #plot(g)
    sis_vpm_simulate(g, 0.4, 0.5, 0.1, 10)
