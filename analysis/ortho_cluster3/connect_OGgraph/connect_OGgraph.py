"""Find connected components of OG graph."""

import os
from src.ortho_cluster.graphs import get_connected_components

# Load graph
graph = {}
with open('../OGs2graph/out/OGgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':')[0] for adj in adjs.split(',')] if adjs else []

# Find connected components
CCs = get_connected_components(graph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/OGconnect.txt', 'w') as outfile:
    for i, CC in enumerate(CCs):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../OGs2graph/OGs2graph.py
    ../OGs2graph/out/OGgraph.tsv
"""