"""Convert OGs to OG graph using cluster overlaps."""

import os
from itertools import combinations

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2gnid[fields['ppid']] = fields['gnid']

# Load OGs
OGid2gnids = {}
with open('../add_paralogs/out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        gnids = {ppid2gnid[node] for edge in fields['edges'].split(',') for node in edge.split(':')}
        OGid2gnids[fields['OGid']] = gnids

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
../add_paralogs/add_paralogs.py
    ../add_paralogs/out/clusters.tsv
"""