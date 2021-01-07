"""Find connected components of ggraph."""

import os
from src.ortho_cluster.DFS import connect

# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = [adj.split(':')[0] for adj in adjs.split(',')]

# Cluster by triangle criterion
CCs = connect(ggraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/gconnect.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../../../src/ortho_cluster/DFS.py
../hits2ggraph/hits2ggraph.py
    ../hits2ggraph/out/ggraph.tsv
"""