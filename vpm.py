# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""
from igraph import *
import random
import numpy as np
import scipy.linalg

dbg = False #Debug flag

try:
    from config import *
except ImportError:
    raise ImportError("Config file unavailable")
    
def get_larg_eig(graph):
    M = len(list(graph.vs))
    l_eig = scipy.linalg.eigh(graph.get_adjacency().data,
                            eigvals_only=True,
                            eigvals=(M-1, M-1))
    return l_eig[0]

def get_eff_strength(graph, b, d, l_eig=None):
    if l_eig==None:
        l_eig = get_larg_eig(graph)
    return l_eig, l_eig*b/d

def immun_random(graph, k):
    N = np.size(graph.vs())
    assert k<=N,'k should be less than N'
    if dbg: print 'initial size: %d'%(N)

    nodes = random.sample(range(N), k)
    graph.delete_vertices(nodes)
    return graph

def immun_highest_degree(graph, k):
    N = np.size(graph.vs())
    assert k<=N,'k should be less than N'

    degrees = [graph.degree(i) for i in range(N)]
    nodes = list(np.argsort(degrees)[-k:])
    graph.delete_vertices(nodes)
    return graph

def immun_highest_degree_iterative(graph, k):
    N = np.size(graph.vs())
    assert k<=N,'k should be less than N'

    for _i in range(k):
        immun_highest_degree(graph, 1)
    return graph

def immun_largest_eigen_vec(graph, k):
    N = np.size(graph.vs())
    assert k<=N,'k should be less than N'

    l_eig = scipy.linalg.eigh(graph.get_adjacency().data,
                            eigvals=(N-1, N-1))
    l_eig_vec =  [i[0] for i in l_eig[1]]
    nodes = list(np.argsort(l_eig_vec)[-k:])
    graph.delete_vertices(nodes)
    return graph

def sis_vpm_simulate(graph, B, D, c, t, immunize=None, k=None,
                     run_simulate=True):
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

    If a policy name is passed as the 'immunize' param, then
    those nodes are immunized at the beginning of the simulation.
    """
    
    immun_dict = {
        'policy_a': lambda: immun_random(graph, k),
        'policy_b': lambda: immun_highest_degree(graph, k),
        'policy_c': lambda: immun_highest_degree_iterative(graph, k),
        'policy_d': lambda: immun_largest_eigen_vec(graph, k)
    }
    if immunize is not None:
        graph = immun_dict[immunize]()

    if not run_simulate:
        larg_eig, eff_strength = get_eff_strength(graph, B1, D1)
        return larg_eig, eff_strength
        
    N = np.size(graph.vs) # Number of nodes
    if dbg: print 'graph size: %d'%(N)
    
    assert 0<=c<=1, ' c should be between 0 and 1'

    infected = set(random.sample(xrange(N),int(c*N)))

    num_infected = [len(infected)]

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



def run_simulation(model, runs, graph, B, D, c, t, immunize=None, k=None):
    """
    Runs the simulation 'runs' # of times and return average nummber
    of infected nodes
    """
    graph_ = graph.copy()
    print 'Running simulation %d times'%(runs)
    if model=='SIS':
        sim_res=range(runs)
        for i in range(runs):
            sim_res[i] = sis_vpm_simulate(graph_, B, D, c, t,
                                          immunize=immunize, k=k)
            graph_ = graph.copy()

        avg_res = []
        for i in range(t):
            avg_res.append(mean([sim_res[j][i] for j in range(runs)]))
        
        print 'Nodes initially infected: %d'%(avg_res[0])
        print 'Avg. Nodes finally infected: %d'%(avg_res[-1])
        if avg_res[0] <= avg_res[-1]:
            print 'Virus has caused an epidemic'
        else:
            print 'Virus epidemic has been prevented'
        return avg_res
        
def num_vaccince_analysis(graph, B, D, immunize, k_list):
    eff_strens=[]
    graph_ = graph.copy()
    for k in k_list:                    
        larg_eig, eff_stren = sis_vpm_simulate(graph_, 
                B, D, None, None, immunize=immunize, 
                k=k, run_simulate=False)
        eff_strens.append(eff_stren)        
        graph_ = graph.copy()
    return eff_strens
    
class Alternating_Networks:
    def system_matrix(self, g1, g2, B, D):
        """
        System-matrix is given by S1xS2, where 
        S1=(1-D)*I + B*A1; 
        S2=(1-D)*I + B*A2;
        (A1 and A2 are the adjacency matrices of the two alternating networks)
        """
        if dbg: print "Calculating S1"
        S1 = (1-D)*np.identity(len(g1.vs))+B*np.array(g1.get_adjacency().data)
        
        if dbg: print "Calculating S2"
        S2 = (1-D)*np.identity(len(g2.vs))+B*np.array(g2.get_adjacency().data)
        
        if dbg: print "Calculating S = S1 X S2"
        return S1.dot(S2)
        
    def get_eff_strength(self, g1, g2, B, D):
        S = self.system_matrix(g1, g2, B, D)
        M = np.shape(S)[0]
        if dbg: print "Calculating largest eigen value"
        l_eig = scipy.linalg.eigh(S, eigvals_only=True, eigvals=(M-1, M-1))
        return l_eig[0]
        
    def sis_vpm_simulate(graphs, B, D, c, t, immunize=None, k=None,
                     run_simulate=True):
        #TODO: Update for graph's'
        """
        immun_dict = {
            'policy_a': lambda: immun_random(graph, k),
            'policy_b': lambda: immun_highest_degree(graph, k),
            'policy_c': lambda: immun_highest_degree_iterative(graph, k),
            'policy_d': lambda: immun_largest_eigen_vec(graph, k)
        }
        if immunize is not None:
            graph = immun_dict[immunize]()
    
        if not run_simulate:
            larg_eig, eff_strength = get_eff_strength(graph, B1, D1)
            return larg_eig, eff_strength
        """
        N = np.size(graphs[0].vs) # Number of nodes
        assert 0<=c<=1, ' c should be between 0 and 1'
        infected = set(random.sample(xrange(N),int(c*N)))
        num_infected = [len(infected)]
    
        #Start simulation
        for _i in range(t):
            graph = graphs[_i%len(graphs)] #Alternate the graphs
            infected_new = set()
            for n in infected:
                nbrs = graph.neighbors(n)
                nbrs_infected = random.sample(nbrs, int(B*len(nbrs)))
                for j in nbrs_infected:
                    infected_new.add(j)
            cured = random.sample(infected, int(math.ceil(D*len(infected))))
    
            for n in infected_new:
                infected.add(n)
    
            for n in cured:
                infected.remove(n)
    
            num_infected.append(len(infected))
        return num_infected
        
    
if __name__=='__main__':
    g=Graph()
    g.add_vertices(10)
    g.add_edges([(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),
                 (4,5),(5,6),(6,7),(7,8),(8,9),(5,9),(5,8)])
    #plot(g)
    print sis_vpm_simulate(g, 0.4, 0.5, 0.1, 10)
