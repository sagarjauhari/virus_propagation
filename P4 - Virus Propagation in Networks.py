# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Virus Propagation in Networks
# Analysis of both **Static Networks** as well as **Alternating Networks** has been described in this report

# <markdowncell>

# #1. Static Networks
# 
# ---

# <codecell>

"""
Created on Sun Nov 24 17:49:17 2013
@author: sagar jauhari
"""
from os.path import join
from igraph import Graph
import scipy.linalg
import matplotlib.pyplot as plt

import vpm
reload(vpm)

pylab.rcParams['figure.figsize'] = 8, 4

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

# <markdowncell>

# ##1. Theoretical Analysis
# For the SIS (susceptible, infected, susceptible) Virus Propagation Model (VPM), with transmission probability β = β 1 , and healing probability δ = δ 1 , calculate the effective strength (s) of the virus on the static contact network provided (static.network). See supplementary material provided for details on the SIS VPM and on how to calculate the effective strength of a virus. Answer the following questions:

# <codecell>

"""
scipy.linalg.eigh
@param eigvals : tuple (lo, hi), optional

    Indexes of the smallest and largest (in ascending 
    order) eigenvalues and corresponding eigenvectors 
    to be returned: 0 <= lo <= hi <= M-1. If omitted, 
    all eigenvalues and eigenvectors are returned.

"""
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

graph = file2igraph(join(DATA_URL,'static.network'))

# <codecell>

larg_eig, eff_strength = get_eff_strength(graph, B1, D1)

print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# ####a. Will the infection spread across the network (i.e., result on an epidemic), or will it die quickly?
# Since the Effective Strength is > 1, the virus will result in an epidemic

# <markdowncell>

# ####b. Keeping δ fixed, analyze how the value of β affects the effective strength of the virus (suggestion: plot your results). What is the minimum transmission probability (β) that results in a network-wide epidemic?

# <codecell>

def plot_change_beta(D, L1):
    b = [i/100.0 for i in range(1,100)]
    s = [get_eff_strength(graph, i, D, l_eig=larg_eig)[1] for i in b]
    plt.plot(b,s,'o-')
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    plt.show()
plot_change_beta(D1, larg_eig)
print 'Minimum Beta that results in network-wide epidemic: %03f' % (D1/larg_eig)

# <markdowncell>

# ####c. Keeping β fixed, analyze how the value of δ affects the effective strength of the virus (suggestion: plot your results). What is the maximum healing probability (δ) that results in a network-wide epidemic?

# <codecell>

def plot_change_delta(B, L1):
    d = [i/100.0 for i in range(1,100)]
    s = [get_eff_strength(graph, B, i, l_eig=larg_eig)[1] for i in d]
    plt.plot(d,s,'o-')
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    plt.show()
    
plot_change_delta(B1, larg_eig)
print 'Maximum Delta that results in network-wide epidemic: %03f' % (B1*larg_eig)

# <markdowncell>

# ####d. Repeat (1), (1a), (1b) and (1c) with β = β 2 , and δ = δ 2 .

# <codecell>

larg_eig, eff_strength = get_eff_strength(graph, B2, D2)
print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# Since effective strength is < 1, there will never be an epidemic

# <codecell>

plot_change_beta(D2, larg_eig)
print 'Minimum Beta that results in network-wide epidemic: %03f' % (D2/larg_eig)

# <codecell>

plot_change_delta(B2, larg_eig)
print 'Maximum Delta that results in network-wide epidemic: %03f' % (B2*larg_eig)

# <markdowncell>

# ##2. Simulation
# Write a program that simulates the propagation of a virus with the SIS VPM, given a static contact
# network, a transmission probability (β), a healing probability (δ), a number of initially infected nodes
# (c), and a number of time steps to run the simulation (t). The initially infected nodes should be chosen
# from a random uniform probability distribution. At each time step, every susceptible (i.e., non-
# infected) node has a β probability of being infected by neighboring infected nodes, and every infected
# node has a δ probability of healing and becoming susceptible again. Your program should also calculate
# the fraction of infected nodes at each time step.
# 
# ####a. Run the simulation program 10 times for the static contact network provided (static.network), with β = β 1 , δ = δ 1 , c = n/10 (n is the number of nodes in the network), and t = 100.

# <codecell>

# simulation program written in file 'vmp.py'
reload(vpm)
result = vpm.run_simulation('SIS', 10, graph, B1, D1, 0.1, 100)

# <markdowncell>

# ####b. Plot the average (over the 10 simulations) fraction of infected nodes at each time step. Did the infection spread across the network, or did it die quickly? Do the results of the simulation agree with your conclusions in (1a)?

# <codecell>

plt.plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# The infection in the network spread and did not die.
# 
# Yes, the results of the simulation agree with the answer of 1(a) - infection resulted in an epidemic!

# <markdowncell>

# ####c. Repeat (2a) and (2b) with β = β 2 , and δ = δ 2 .

# <codecell>

result = vpm.run_simulation('SIS', 10, graph, B2, D2, 0.1, 100)
plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# The infection in the network quickly died out.
# 
# Yes, the results of the simulation agree with the answer of 1(d) - infection did not result in an epidemic but quickly died out.

# <markdowncell>

# ##3. Immunization Policy
# Write a program that implements an immunization policy to prevent the virus from spreading across
# the network. Given a number of available vaccines (k) and a contact network, your program should
# select k nodes to immunize. The immunized nodes (and their incident edges) are then removed from
# the contact network.  

# <markdowncell>

# ####a. What do you think would be the optimal immunization policy? What would be its time complexity? Would it be reasonable to implement this policy? Justify.
# - The optimal immunization policy would be to select k nodes whose removal will cause the largest drop in the value of λ 1 for immunization (largest eigen-drop). But the time complexity of this algorithm is very high because we'll have to calculate the eigen drop for each of the possible subsets of k nodes.
# - From amongst Policy A, B, C, D, Policy D should be used because it caused the largest drop in the eigen value (new eigen value after immunization was 0.583) and prevented the epidemic. The compuational complexity is reasonable since not highly intensive calculations are required for the strategy.

# <markdowncell>

# ####For your program, use the following heuristic immunization policies:
# 
# - Policy A: Select k random nodes for immunization.
# - Policy B: Select the k nodes with highest degree for immunization.
# - Policy C: Select the node with the highest degree for immunization. Remove this node (and its incident edges) from the contact network. Repeat until all vaccines are administered.
# - Policy D: Find the eigenvector corresponding to the largest eigenvalue of the contact network’s adjacency matrix. Find the k largest (absolute) values in the eigenvector. Select the k nodes at the corresponding positions in the eigenvector.
# 
# ####For each heuristic immunization policy (A, B, C, and D) and for the static contact network provided
# (static.network), answer the following questions:
# 
# b. What do you think is the intuition behind this heuristic?
# 
# c. Write a pseudocode for this heuristic immunization policy. What is its time complexity?
# 
# d. Given k = k 1 , β = β 1 , and δ = δ 1 , calculate the effective strength (s) of the virus on the immunized
# contact network (i.e., contact network without immunized nodes). Did the immunization policy
# prevented a network-wide epidemic?
# 
# e. Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.
# 
# f. Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <markdowncell>

# ### Policy A
# 
# ----
# #### Intuition
# Select k nodes randomly. If we use this strategy on a lot of graphs, at least some times we will be able to contain the epidemic  

# <markdowncell>

# #### Pseudo Code
# <pre>
# def immun_random(graph, k):
#     N = Number of nodes in graph
#     nodes = random.sample(range(N), k)
#     graph.delete_vertices(nodes)
#     return graph
# </pre>
# #### Effective Strength

# <codecell>

reload(vpm)
larg_eig, eff_strength = vpm.sis_vpm_simulate(graph, 
                B1, D1, 
                None, None,
                immunize='policy_a', 
                k=K1,
                run_simulate=False)
print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus is still very high and the epidemic will not be contained.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

k_list = [4500 + i*100 for i in range(10)]
strens = vpm.num_vaccince_analysis(graph, B1, D1,'policy_a',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that a very large number of vaccines will be required to prevent the epidemic: ~5200 - 5300

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

reload(vpm)
result_a = vpm.run_simulation('SIS', 10, graph,
                            B1, D1, 0.1, 100, 
                            immunize='policy_a', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy B
# 
# ----
# #### Intuition
# The idea is to remove nodes of highest degree. By doing this, we're making sure that the average degree of the nodes of the graph decreases so we make it difficult for the virus to propagate

# <markdowncell>

# #### Pseudo Code
# <pre>
# def immun_highest_degree(graph, k):
#     N = Number of nodes in graph
#     degrees = [graph.degree(i) for i in range(N)]
#     nodes = list(np.argsort(degrees)[-k:])
#     graph.delete_vertices(nodes)
#     return graph
# </pre>
# #### Effective Strength

# <codecell>

blarg_eig, eff_strength = vpm.sis_vpm_simulate(graph, 
                B1, D1, 
                None, None,
                immunize='policy_b', 
                k=K1,
                run_simulate=False)
print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization (~1) and the virus should be contained.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

k_list = [i*50 for i in range(4,8)]
strens = vpm.num_vaccince_analysis(graph, B1, D1,'policy_b',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~200 - 250

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

result = vpm.run_simulation('SIS', 10, graph,
                            B1, D1, 0.1, 100, 
                            immunize='policy_b', k=K1)
plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy C
# 
# ---
# #### Intuition
# Unlike the previous method where 'k' highest degrees are removed all together, this policy follows an iterative approach. After removing the node of highest degree from the graph, it recalculates the degree of each node to find out the new node of highest degree. This is computationally more expensive than the previous policy by theoretically better in reducing the average connectivity of the graph 

# <markdowncell>

# #### Pseudo Code
# <pre>
# def immun_highest_degree_iterative(graph, k):
#     N = Number of nodes in graph
#     for _i in range(k):
#         Remove highest degree node
#     return graph
# </pre>
# #### Effective Strength

# <codecell>

clarg_eig, eff_strength = vpm.sis_vpm_simulate(graph, 
                B1, D1, 
                None, None,
                immunize='policy_c', 
                k=K1,
                run_simulate=False)
print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization ( < 1 ) and the epidemic would be prevented.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

k_list = [i*10 for i in range(14,25)]
strens = vpm.num_vaccince_analysis(graph, B1, D1,'policy_c',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~220 - 230

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

result = vpm.run_simulation('SIS', 10, graph,
                            B1, D1, 0.1, 100, 
                            immunize='policy_c', k=K1)
plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy D
# 
# ---
# #### Intuition
# This policy is a very simplified version of the NetShield algorithm. The idea is to remove then nodes corresponding to the largest values in the eigen vector of the largest eigen value. By removing the nodes, it is expected that the largest eigen value and thus decrease in the effective strength of the virus

# <markdowncell>

# #### Pseudo Code
# <pre>
# def immun_largest_eigen_vec(graph, k):
#     N = Number of nodes in graph
#     l_eig = Largest igen value
#     l_eig_vec =  Eigen vector of largest eigen value
#     nodes = Find index of 'k' largest values in l_eig_vec
#     graph.delete_vertices(nodes)
#     return graph
# </pre>
# #### Effective Strength

# <codecell>

larg_eig, eff_strength = vpm.sis_vpm_simulate(graph, 
                B1, D1, 
                None, None,
                immunize='policy_d', 
                k=K1,
                run_simulate=False)
print 'Largest Eigen Value: %0.3f'% (larg_eig)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization ( < 1 ) and the epidemic would be prevented.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

k_list = [i*10 for i in range(35, 44)]
strens = vpm.num_vaccince_analysis(graph, B1, D1,'policy_d',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~420 - 430

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

result_a = vpm.run_simulation('SIS', 10, graph,
                            B1, D1, 0.1, 100, 
                            immunize='policy_d', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# #2. Virus Propagation on Time-Varying Networks
# 
# ---

# <markdowncell>

# ##1. Theoretical Analysis
# 
# Repeat problem (1) using the set of 2 alternating contact networks provided (alternating1.network and
# alternating2.network). See supplementary material provided for details on how to calculate the
# effective strength of a virus for time-varying contact networks. Answer questions (1), (1a), (1b), (1c),
# and (1d).
# 
# For the SIS (susceptible, infected, susceptible) Virus Propagation Model (VPM), with transmission probability β = β 1 , and healing probability δ = δ 1 , calculate the effective strength (s) of the virus on the static contact network provided (static.network). See supplementary material provided for details on the SIS VPM and on how to calculate the effective strength of a virus. Answer the following questions:

# <codecell>

graph1 = file2igraph(join(DATA_URL,'alternating1.network'))
graph2 = file2igraph(join(DATA_URL,'alternating1.network'))
graphs=[graph1, graph2]
reload(vpm)
acn = vpm.Alternating_Networks()

# <codecell>

eff_strength = acn.get_eff_strength(graph1, graph2, B1, D1)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# ####a. Will the infection spread across the network (i.e., result on an epidemic), or will it die quickly?
# Since the Effective Strength is > 1, the virus will result in an epidemic

# <markdowncell>

# ####b. Keeping δ fixed, analyze how the value of β affects the effective strength of the virus (suggestion: plot your results). What is the minimum transmission probability (β) that results in a network-wide epidemic?

# <codecell>

def plot_change_beta(g1, g2, D, b):
    s = [acn.get_eff_strength(g1, g2, i, D) for i in b]
    plt.plot(b,s,'o-')
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    plt.show()
plot_change_beta(graph1, graph2, D1, 
    [i/100.0 for i in range(1,10)])

# <markdowncell>

# The minimum value of B required for propagation is very less: ~0.15

# <markdowncell>

# ####c. Keeping β fixed, analyze how the value of δ affects the effective strength of the virus (suggestion: plot your results). What is the maximum healing probability (δ) that results in a network-wide epidemic?

# <codecell>

def plot_change_delta(g1, g2, B):
    d = [i/1000.0 for i in range(995,1000)]
    s = [acn.get_eff_strength(g1, g2, B, i) for i in d]
    plt.plot(d, s,'o-')
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    plt.show()
    return s
plot_change_delta(graph1, graph2, B1)

# <markdowncell>

# Plot shows that the maximum healing probability which will cause an epidemic is very high: > 0.999

# <markdowncell>

# ####d. Repeat (1), (1a), (1b) and (1c) with β = β 2 , and δ = δ 2 .

# <codecell>

eff_strength = acn.get_eff_strength(graph1, graph2, B2, D2)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength is < 1. So the virus should *NOT* cause an epidemic!

# <codecell>

def plot_change_beta(g1, g2, D, b):
    s = [acn.get_eff_strength(g1, g2, i, D) for i in b]
    plt.plot(b,s,'o-')
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    plt.show()
plot_change_beta(graph1, graph2, D2)

# <markdowncell>


# <codecell>

def plot_change_delta(g1, g2, B):
    d = [i/1000.0 for i in range(995,1000)]
    s = [acn.get_eff_strength(g1, g2, B, i) for i in d]
    plt.plot(d, s,'o-')
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    plt.show()
    return s
plot_change_delta(graph1, graph2, B2)

# <markdowncell>


# <markdowncell>

# ##2. Simulation
# Extend your virus propagation simulation program from problem (2) to allow time-varying graphs. Your
# modified program will be given a set of T alternating contact networks, a transmission probability (β), a
# healing probability (δ), a number of initially infected nodes (c), and a number of time steps to run the
# simulation (t). Your program should alternate between contact networks at each time step. That is, at a
# given time step t i , for i ∈ [0, t), your program should simulate the propagation of the virus on
# alternating contact network (t i mod T) + 1.
# For the modified virus propagation simulation program, and the set of 2 alternating contact networks
# provided (alternating1.network and alternating2.network), answer questions (2a), (2b), and (2c).
# 
# ####a. Run the simulation program 10 times for the static contact network provided (static.network), with β = β 1 , δ = δ 1 , c = n/10 (n is the number of nodes in the network), and t = 100.

# <codecell>

# simulation program written in file 'vmp.py'
reload(vpm)
acn = vpm.Alternating_Networks()
graphs = [graph1, graph2]
result = acn.run_simulation('SIS', 10, graphs, B1, D1, 0.1, 100)

# <markdowncell>

# ####b. Plot the average (over the 10 simulations) fraction of infected nodes at each time step. Did the infection spread across the network, or did it die quickly? Do the results of the simulation agree with your conclusions in (1a)?

# <codecell>

plt.plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the results agree with conclusion in 1(a)

# <markdowncell>

# ####c. Repeat (2a) and (2b) with β = β 2 , and δ = δ 2 .

# <codecell>

result = acn.run_simulation('SIS', 10, graphs, B2, D2, 0.1, 100)
plot(result)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# The infection in the network quickly died out.
# 
# Yes, the results of the simulation agree with the answer of 1(d) - infection did not result in an epidemic but quickly died out.

# <markdowncell>

# ## 3. Immunization
# 
# 1. Extend your immunization policy program from problem (3) to allow time-varying graphs. Your
# modified program will be given a set of T alternating contact networks and a number of available
# vaccines (k). Answer question (3a).
# For your program, use the following heuristic immunization policies:
# - Policy A: Select k random nodes for immunization.
# - Policy B: Select the k nodes with highest average degree across all alternating contact networks for
# immunization.
# - Policy C: Select the node with the highest degree out of all the alternating contact networks for
# immunization. Remove this node (and its incident edges) from all the alternating contact networks.
# Repeat until all vaccines are administered.
# - Policy D: Find the eigenvector corresponding to the largest eigenvalue of the system-matrix (see
# supplementary materials for details on the system-matrix). Find the k largest (absolute) values in
# the eigenvector. Select the k nodes at the corresponding positions in the eigenvector.For each heuristic immunization policy (A, B, C, and D) and for the set of 2 alternating contact networks
# provided (alternating1.network and alternating2.network), answer questions (3b), (3c), (3d), and (3e).
# 
# For each heuristic immunization policy (A, B, C, and D) and for the set of 2 alternating contact networks
# provided (alternating1.network and alternating2.network), answer questions (3b), (3c), (3d), and (3e).

# <markdowncell>

# ####a. What do you think would be the optimal immunization policy? What would be its time complexity? Would it be reasonable to implement this policy? Justify.
# - The optimal immunization policy would be to select k nodes whose removal will cause the largest drop in the value of λ 1 for immunization (largest eigen-drop of the System Matrix). But the time complexity of this algorithm is very high because we'll have to calculate the eigen drop for each of the possible subsets of k nodes.
# - From amongst Policy A, B, C, D, Policy D should be used because it caused the largest drop in the eigen value (new eigen value after immunization was 0.583) and prevented the epidemic. The compuational complexity is reasonable since not highly intensive calculations are required for the strategy.

# <markdowncell>

# ### Policy A
# 
# ----
# #### Intuition
# Select k nodes randomly. If we use this strategy on a lot of graphs, at least some times we will be able to contain the epidemic  

# <markdowncell>

# #### Pseudo Code
# <pre>
#     def immun_random(self, graphs, k):
#         N = np.size(graphs[0].vs())
#         nodes = random.sample(range(N), k)
#         for g in graphs:
#             g.delete_vertices(nodes)
#         return graphs
# </pre>
# #### Effective Strength

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
eff_strength = acn.sis_vpm_simulate(graphs, 
                B1, D1, 
                None, None,
                immunize='policy_a', 
                k=K1,
                run_simulate=False)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus is still very high and the epidemic will not be contained.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()

k_list = [4500 + i*200 for i in range(10, 16)]
strens = acn.num_vaccince_analysis(graphs, B1, D1,'policy_a',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~7500

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

graphs=[graph1, graph2]
[len(i.vs) for i in graphs]

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
result_a = acn.run_simulation('SIS', 10, graphs,
                            B1, D1, 0.1, 100, 
                            immunize='policy_a', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy B
# 
# ----
# #### Intuition
# The idea is to remove nodes of highest degree. By doing this, we're making sure that the average degree of the nodes of the graph decreases so we make it difficult for the virus to propagate

# <markdowncell>

# #### Pseudo Code
# <pre>
#     def immun_highest_degree(self, graphs, k):
#         degrees = [mean(graphs[0].degree(i), 
#                         graphs[1].degree(1)) for i in range(N)]
#         nodes = list(np.argsort(degrees)[-k:])
#         for g in graphs:        
#             g.delete_vertices(nodes)
#         return graphs
# </pre>
# #### Effective Strength

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
eff_strength = acn.sis_vpm_simulate(graphs, 
                B1, D1, 
                None, None,
                immunize='policy_b', 
                k=K1,
                run_simulate=False)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization (~1) and the virus should be contained.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()

k_list = [i*50 for i in range(4, 7)]
strens = acn.num_vaccince_analysis(graphs, B1, D1,'policy_b',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~260

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
result_a = acn.run_simulation('SIS', 10, graphs,
                            B1, D1, 0.1, 100, 
                            immunize='policy_b', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy C
# 
# ---
# #### Intuition
# Unlike the previous method where 'k' highest degrees are removed all together, this policy follows an iterative approach. After removing the node of highest degree from the graph, it recalculates the degree of each node to find out the new node of highest degree. This is computationally more expensive than the previous policy by theoretically better in reducing the average connectivity of the graph 

# <markdowncell>

# #### Pseudo Code
# <pre>
#     def immun_highest_degree_iterative(self, graphs, k):
#         for _i in range(k):
#             immun_highest_degree(graphs, 1)
#         return graphs
# </pre>
# #### Effective Strength

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
eff_strength = acn.sis_vpm_simulate(graphs, 
                B1, D1, 
                None, None,
                immunize='policy_c', 
                k=K1,
                run_simulate=False)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization ( < 1 ) and the epidemic would be prevented.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()

k_list = [i*100 for i in range(1, 4)]
strens = acn.num_vaccince_analysis(graphs, B1, D1,'policy_c',k_list)
plot(k_list, strens, 'o-')
plt.xlabel('Immunized Nodes')
plt.ylabel('Effective strength of virus')

# <markdowncell>

# Analysis shows that the number of vaccines required to prevent the epidemic: ~250

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
result_a = acn.run_simulation('SIS', 10, graphs,
                            B1, D1, 0.1, 100, 
                            immunize='policy_c', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

# <markdowncell>

# ### Policy D
# 
# ---
# #### Intuition
# This policy is a very simplified version of the NetShield algorithm. The idea is to remove then nodes corresponding to the largest values in the eigen vector of the largest eigen value. By removing the nodes, it is expected that the largest eigen value and thus decrease in the effective strength of the virus

# <markdowncell>

# #### Pseudo Code
# <pre>
#     def immun_largest_eigen_vec(self, graphs, k, B, D):
#         S = system_matrix(self, graphs[0], graphs[1], B, D)
#         M=np.shape(S)[0]
#     
#         l_eig = scipy.linalg.eigh(S, eigvals=(M-1, M-1))
#         l_eig_vec =  [i[0] for i in l_eig[1]]
#         nodes = list(np.argsort(l_eig_vec)[-k:])
#         for g in graphs:        
#             g.delete_vertices(nodes)
#         return graphs
# </pre>
# #### Effective Strength

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
eff_strength = acn.sis_vpm_simulate(graphs, 
                B1, D1, 
                None, None,
                immunize='policy_d', 
                k=K1,
                run_simulate=False)
print 'Effective Strength: %0.3f'%(eff_strength)

# <markdowncell>

# The effective strength of the virus has reduced drastically after immunization ( < 1 ) and the epidemic would be prevented.

# <markdowncell>

# #### Analysis of k
# Keeping β and δ fixed, analyze how the value of k affects the effective strength of the virus on
# the immunized contact network (suggestion: plot your results). Estimate the minimum number
# of vaccines necessary to prevent a network-wide epidemic.

# <markdowncell>

# **This is taking more than an hour to calculate!**

# <markdowncell>

# #### Simulation of immunized network
# Given k = k 1 , β = β 1 , δ = δ 1 , c = n/10, and t = 100, run the simulation from problem (2) for the
# immunized contact network 10 times. Plot the average fraction of infected nodes at each time
# step. Do the results of the simulation agree with your conclusions in (3d)?

# <codecell>

reload(vpm)
acn = vpm.Alternating_Networks()
result_a = acn.run_simulation('SIS', 10, graphs,
                            B1, D1, 0.1, 100, 
                            immunize='policy_d', k=K1)
plot(result_a)
plt.xlabel('Time')
plt.ylabel('Infected Nodes')

# <markdowncell>

# Yes, the simulataion results agree with the conclusion based on effective strength of virus

