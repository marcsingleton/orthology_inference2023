"""Cluster ggraph on triangle criterion."""

import json
import os
from gtriDFS import cluster

# Load best hits graph
with open('../blast2ggraph/out/ggraph.json') as file:
    ggraph = json.load(file)

# Remove non-reciprocal hits
num_del = 0
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
            num_del += 1

    # Remove non-reciprocal hits after initial loop is completed to not modify list during loop
    for del_key in del_keys:
        del adj_gns[del_key]

# Cluster by triangle criterion
OGs = cluster(ggraph)

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Write clusters to file
with open('out/gclusters.json', 'w') as outfile:
    json.dump(OGs, outfile, indent=1)

"""
DEPENDENCIES
../blast2ggraph/blast2ggraph.py
    ../blast2ggraph/out/ggraph.json
./gtriDFS.py
"""