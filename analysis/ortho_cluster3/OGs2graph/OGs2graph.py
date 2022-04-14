"""Convert OGs to OG graph using cluster overlaps."""

import os
from itertools import combinations

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.rstrip('\n').split('\t')
        ppid2gnid[ppid] = gnid

# Load OGs
OGid2gnids = {}
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, OGid, _, edges = line.rstrip('\n').split('\t')
        gnids = {ppid2gnid[node] for edge in edges.split(',') for node in edge.split(':')}
        OGid2gnids[OGid] = gnids

# Make OG graph
graph = {OGid: set() for OGid in OGid2gnids}
for OGid1, OGid2 in combinations(OGid2gnids, 2):
    gnids1, gnids2 = OGid2gnids[OGid1], OGid2gnids[OGid2]
    if len(gnids1 & gnids2) / min(len(gnids1), len(gnids2)) >= 0.5:
        graph[OGid1].add(OGid2)
        graph[OGid2].add(OGid1)

# Write OG graph to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/OG_graph.tsv', 'w') as file:
    for node, adjs in graph.items():
        file.write(node + '\t' + ','.join(adjs) + '\n')

"""
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../cluster4+_graph/cluster.py
    ../cluster4+_graph/out/4clique/clusters.tsv
"""