"""Cluster pgraph on triangle criterion."""

import os
from src.ortho_cluster.triDFS import cluster

# Load pgraph
pgraph = {}
with open('../hsps2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = set(adjs.split(','))

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

# Cluster by triangle criterion
OGs = cluster(rpgraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/pclusters.txt', 'w') as outfile:
    for i, OG in enumerate(OGs):
        OGid = hex(i)[2:].zfill(4)
        outfile.write(OGid + ':' + '\t'.join([f'{node1},{node2}' for node1, node2 in OG]) + '\n')

"""
DEPENDENCIES
../hsps2pgraph/hsps2pgraph.py
    ../hsps2pgraph/out/pgraph.tsv
"""