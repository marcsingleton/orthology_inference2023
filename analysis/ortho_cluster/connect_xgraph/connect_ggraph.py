"""Find connected components of ggraph."""

import os
import json
from DFS import connect

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

# Remove substructure
for node, adj_gns in ggraph.items():
    ggraph[node] = list(adj_gns)

# Cluster by triangle criterion
CCs = connect(ggraph)

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Write clusters to file
with open('out/gconnect.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../blast2pgraph/blast2ggraph.py
    ../blast2pgraph/out/ggraph.json
./DFS.py
"""