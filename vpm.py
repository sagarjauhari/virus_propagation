# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:49:17 2013

@author: sagar jauhari
"""

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

    The initially infected nodes should be chosen from a random uniform
    probability distribution. At each time step, every susceptible (i.e., non-
    infected) node has a β probability of being infected by neighboring
    infected nodes, and every infected node has a δ probability of healing and
    becoming susceptible again. The program also calculates the fraction
    of infected nodes at each time step.
    """
    pass
