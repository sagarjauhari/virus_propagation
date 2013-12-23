Virus Propagation
=================

Virus Propagation on Static Contact Networks and Alternating Contact Networks

View [IPython Notebook](http://nbviewer.ipython.org/github/sagarjauhari/virus_propagation/blob/master/P4%20-%20Virus%20Propagation%20in%20Networks.ipynb)

- Epidemic Threshold: Theoretical analysis of virus propagation based on effective strength [1]
- Simulation of virus propagation across the network
- Analylsis of different immunization policies
    - Policy A: Select k random nodes for immunization.
    - Policy B: Select the k nodes with highest degree for immunization.
    - Policy C: Select the node with the highest degree for immunization. Remove this node (and its incident edges) from the contact network. Repeat until all vaccines are administered.
    - Policy D: Find the eigenvector corresponding to the largest eigenvalue of the contact network’s adjacency matrix. Find the k largest (absolute) values in the eigenvector. Select the k nodes at the corresponding positions in the eigenvector.



##References
[1] Prakash, B. Aditya, et al. "Got the Flu (or Mumps)? Check the Eigenvalue!." arXiv preprint arXiv:1004.0060 (2010).
