# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""
from igraph import *
import random
import numpy as np
import scipy.linalg

try:
    from config import *
except ImportError:
    raise ImportError("Config file unavailable")

def sis_vpm_simulate(graph, B, D, c, t, immunize=None):
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

    If there is a list of nodes passed as the 'immunize' param, then
    those nodes are immunized at the beginning of the simulation.
    """
    assert 0<=c<=1, ' c should be between 0 and 1'

    if immunize is not None:
        graph.delete_vertices(immunize)

    N = len(list(graph.vs)) # Number of nodes
    infected = set(random.sample(xrange(N),int(c*N)))

    num_infected = []

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

        num_infected.append(len(infected))
    return num_infected

def immun_random(graph, k):
    N = size(graph.vs())
    assert k<=N,'k should be less than N'

    nodes = random.sample(range(N), k)
    graph.delete_vertices(nodes)
    return graph

def immun_highest_degree(graph, k):
    N = size(graph.vs())
    assert k<=N,'k should be less than N'

    degrees = [graph.degree(i) for i in range(N)]
    nodes = list(np.argsort(degrees)[-k:])
    graph.delete_vertices(nodes)
    return graph

def immun_highest_degree_iterative(graph, k):
    N = size(graph.vs())
    assert k<=N,'k should be less than N'

    for _i in range(k):
        node = immun_highest_degree(graph, 1)
        gg.delete_vertices(node)

    return graph


def immun_largest_eigen_vec(graph, k):
    N = size(graph.vs())
    assert k<=N,'k should be less than N'

    l_eig = scipy.linalg.eigh(graph.get_adjacency().data,
                            eigvals=(N-1, N-1))
    l_eig_vec =  [i[0] for i in l_eig[1]]
    nodes = list(np.argsort(l_eig_vec)[-k:])
    graph.delete_vertices(nodes)
    return graph

def run_simulation(model, runs, graph, B, D, c, t):
    """
    Runs the simulation 'runs' # of times and return average nummber
    of infected nodes
    """
    if model=='SIS':
        sim_res=range(runs)
        for i in range(runs):
            sim_res[i] = sis_vpm_simulate(graph, B, D, c, t)

        avg_res = []
        for i in range(t):
            avg_res.append(mean([sim_res[j][i] for j in range(runs)]))
        return avg_res

if __name__=='__main__':
    g=Graph()
    g.add_vertices(10)
    g.add_edges([(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),
                 (4,5),(5,6),(6,7),(7,8),(8,9),(5,9),(5,8)])
    #plot(g)
    print sis_vpm_simulate(g, 0.4, 0.5, 0.1, 10)
