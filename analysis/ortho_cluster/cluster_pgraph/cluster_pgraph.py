"""Cluster pgraph on triangle criterion."""

import os
from ptriDFS import cluster

# Parse best hits as graph
pgraph = {}
with open('../blast2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Remove non-reciprocal hits
for node, adjs in pgraph.items():
    # Search current node for non-reciprocal hits
    adj_idxs = []
    for adj_idx, adj in enumerate(adjs):
        try:  # Cannot test with "in" easily since it assumes the node is in the graph in the first place
            pgraph[adj].index(node)
        except (KeyError, ValueError):  # KeyError from adj not in tgraph; ValueError from node not in adjs
            adj_idxs.append(adj_idx)

    # Remove non-reciprocal hits after initial loop is completed to not modify list during loop
    for offset, adj_idx in enumerate(adj_idxs):
        del adjs[adj_idx - offset]

# Cluster by triangle criterion
OGs = cluster(pgraph)

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Write clusters to file
with open('out/pclusters.tsv', 'w') as outfile:
    for OG in OGs:
        outfile.write(','.join(OG) + '\n')

"""
DEPENDENCIES
../blast2pgraph/blast2pgraph.py
    ../blast2pgraph/out/pgraph.tsv
./ptriDFS.py
"""