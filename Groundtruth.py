
import networkx as nx
import numpy as np
import sys

lines = []
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

nodes = []
edges = []

count = 0

length = len(lines)
for line in lines:
    count += 1
    print(line)
    if count == 2:
        nodes = line.strip().split(';')
        print(nodes)
    if count > 4 and count < (length):
        temp = []
        temp = line.strip().split(' ')
        
        edges.append((temp[1], temp[3]))
        print(edges)

G = nx.DiGraph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)

import matplotlib.pyplot as plt
nx.draw(G, with_labels=True)
plt.show()

# https://fentechsolutions.github.io/CausalDiscoveryToolbox/html/metrics.html

A = nx.to_numpy_matrix(G)

filename='adjacency_matrix.txt'
np.savetxt(filename, A, fmt='%d')

