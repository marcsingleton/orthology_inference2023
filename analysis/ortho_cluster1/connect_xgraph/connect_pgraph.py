"""Find connected components of pgraph."""

import os
from src.ortho_cluster.DFS import connect

# Load pgraph
pgraph = {}
with open('../hsps2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Make reciprocal
rpgraph = {}
for qppid, sppids in pgraph.items():
    for sppid in sppids:
        try:
            r = qppid in pgraph[sppid]
        except KeyError:
            r = False

        if r:
            try:
                rpgraph[qppid].add(sppid)
            except KeyError:
                rpgraph[qppid] = {sppid}

# Find connected components
CCs = connect(rpgraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/pconnect.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../hsps2pgraph/hsps2pgraph.py
    ../hsps2pgraph/out/pgraph.tsv
"""