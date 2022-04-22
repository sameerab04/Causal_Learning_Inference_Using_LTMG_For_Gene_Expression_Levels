
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

G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)

gene_names = []
f = open(sys.argv[2])
for row in f:
    gene_names.append(row.strip())

import matplotlib.pyplot as plt
nx.draw(G, with_labels=True, pos = nx.circular_layout(G))

# https://fentechsolutions.github.io/CausalDiscoveryToolbox/html/metrics.html

A = nx.to_numpy_matrix(G, gene_names)

print(G.nodes)
filename ='adjacency_matrix2.txt'

np.savetxt(filename, A, fmt='%d')

plt.show()
