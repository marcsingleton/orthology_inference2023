"""Cluster ggraph on triangle criterion."""

import os
from triDFS import cluster

# Parse best hits as graph
ggraph = {}
with open('../blast2ggraph/out/ggraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        if node != 'null':  # Remove None first to prevent recognition later
            ggraph[node] = adjs.split(',')

# Remove non-reciprocal hits
for node, adjs in ggraph.items():
    # Search current node for non-reciprocal hits
    adj_idxs = []
    for adj_idx, adj in enumerate(adjs):
        try:  # Cannot test with "in" easily since it assumes the node is in the graph in the first place
            ggraph[adj].index(node)
        except (KeyError, ValueError):  # KeyError from adj not in pgraph; ValueError from node not in adjs
            adj_idxs.append(adj_idx)

    # Remove non-reciprocal hits after initial loop is completed to not modify list during loop
    for offset, adj_idx in enumerate(adj_idxs):
        del adjs[adj_idx - offset]

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
../blast2ggraph/blast2ggraph.py
    ../blast2ggraph/out/ggraph.json
./triDFS.py
"""