
import networkx as nx
import numpy as np
import sys

lines = []
f = open(sys.argv[1])
for row in f:
    lines.append((row))

nodes = []
edges = []

count = 1

length = len(lines)
for line in lines:
    if count >= 2:
        nodes = line.strip().split(',')
        if nodes[0] != nodes [1]:
            edges.append((nodes[0],nodes[1]))
    count += 1

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

