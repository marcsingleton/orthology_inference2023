"""Cluster ggraph on triangle criterion."""

import json
import os
from triDFS import cluster

# Load best hits graph
with open('../blast2ggraph/out/ggraph.json') as file:
    ggraph = json.load(file)

# Remove non-reciprocal hits
for node, adj_gns in ggraph.items():
    # Remove None first to prevent recognition
    if 'null' in adj_gns:
        del adj_gns['null']

    # Search current node for non-reciprocal hits
    del_keys = []
    for adj_gn in adj_gns:
        try:
            ggraph[adj_gn][node]
        except KeyError:
            del_keys.append(adj_gn)

    # Remove non-reciprocal hits after initial loop is completed to not modify list during loop
    for del_key in del_keys:
        del adj_gns[del_key]

# Cluster by triangle criterion
OGs = cluster({node: list(adj_gns) for node, adj_gns in ggraph.items()})

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