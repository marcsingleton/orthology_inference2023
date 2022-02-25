"""Convert OGs to OG graph using cluster overlaps."""

import os
from itertools import combinations


# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load OGs
OGid2gnids = {}
with open('../clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = {ppid2gnid[node] for edge in edges.split('\t') for node in edge.split(',')}
        OGid2gnids[OGid] = gnids

# Make OGgraph
graph = {OGid: set() for OGid in OGid2gnids}
for OGid1, OGid2 in combinations(OGid2gnids, 2):
    gnids1, gnids2 = OGid2gnids[OGid1], OGid2gnids[OGid2]
    if len(gnids1 & gnids2) / min(len(gnids1), len(gnids2)) >= 0.5:
        graph[OGid1].add(OGid2)
        graph[OGid2].add(OGid1)

# Write OGgraph to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/OGgraph.tsv', 'w') as file:
    for node, adjs in graph.items():
        file.write(node + '\t' + ','.join(adjs) + '\n')

"""
DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_pcommunity/clique4+_pcommunity2.py
    ../clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
"""