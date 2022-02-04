"""Find connected components of pgraph."""

import os
from src.ortho_cluster.DFS import connect

# Load graph
graph = {}
with open('../hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':')[0] for adj in adjs.split(',')]

# Find connected components
CCs = connect(graph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/pconnect2.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../hits2pgraph/hits2pgraph2.py
    ../hits2pgraph/out/pgraph2.tsv
"""