"""Cluster ggraph on triangle criterion."""

import os
from triDFS import cluster

# Parse best hits as graph
ggraph = {}
with open('../make_xreciprocal/out/ggraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = adjs.split(',')

# Cluster by triangle criterion
OGs = cluster(ggraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/gclusters.txt', 'w') as outfile:
    for i, OG in enumerate(OGs):
        OGid = hex(i)[2:].zfill(4)
        outfile.write(OGid + ':' + '\t'.join([f'{node1},{node2}' for node1, node2 in OG]) + '\n')

"""
DEPENDENCIES
../make_xreciprocal/make_greciprocal.py
    ../make_xreciprocal/out/ggraph.tsv
./triDFS.py
"""