"""Find connected components of ggraph."""

import os
from DFS import connect

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
    ../blast2pgraph/out/ggraph.tsv
./DFS.py
"""