"""Filter OGs to representatives by uniqueness criteria and species number."""

import os

import matplotlib.pyplot as plt
import pandas as pd


def spid_filter(spids):
    conditions = [({'dnov', 'dvir'}, 1),
                  ({'dmoj', 'dnav'}, 1),
                  ({'dinn', 'dgri', 'dhyd'}, 2),
                  ({'dgua', 'dsob'}, 1),
                  ({'dbip', 'dana'}, 1),
                  ({'dser', 'dkik'}, 1),
                  ({'dele', 'dfik'}, 1),
                  ({'dtak', 'dbia'}, 1),
                  ({'dsuz', 'dspu'}, 1),
                  ({'dsan', 'dyak'}, 1),
                  ({'dmel'}, 1),
                  ({'dmau', 'dsim', 'dsec'}, 1)]
    return all([len(spids & group) >= num for group, num in conditions])


spid_min = 20

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2data[ppid] = (gnid, spid, sqid)

# Load graph
graph = {}
with open('../../ortho_cluster3/hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        graph[node] = bitscores

# Load GGs
OGid2GGid = {}
with open('../../ortho_cluster3/connect_OG_graph/out/components.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        GGid, OGids = line.rstrip().split('\t')
        for OGid in OGids.split(','):
            OGid2GGid[OGid] = GGid

# Load OGs
GGid2OGs = {}
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip().split('\t')
        edges = [edge.split(':') for edge in edges.split(',')]
        ppids = {node for edge in edges for node in edge}
        gnids = {ppid2data[ppid][0] for ppid in ppids}
        spids = {ppid2data[ppid][1] for ppid in ppids}
        sqids = {ppid2data[ppid][2] for ppid in ppids}
        bitscore = 0
        for node1, node2 in edges:
            bitscore += graph[node1][node2] + graph[node2][node1]

        if len(gnids) == len(spids) >= spid_min and spid_filter(spids):
            GGid = OGid2GGid[OGid]
            record = {'component_id': component_id, 'OGid': OGid, 'GGid': GGid,
                      'bitscore': round(bitscore, 1), 'edgenum': len(edges),
                      'ppidnum': len(ppids), 'sqidnum': len(sqids),
                      'gnidnum': len(gnids), 'spidnum': len(spids)}
            try:
                GGid2OGs[GGid].append(record)
            except KeyError:
                GGid2OGs[GGid] = [record]

rows = []
for OGs in GGid2OGs.values():
    OG = max(OGs, key=lambda x: (x['spidnum'], x['bitscore']))  # Most species, highest bitscore
    rows.append(OG)
OGs = pd.DataFrame(rows)

if not os.path.exists('out/'):
    os.mkdir('out/')

print('Total filtered OGs:', len(OGs))
counts = OGs['gnidnum'].value_counts()
plt.bar(counts.index, counts.values)
plt.xlabel('Number of genes in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/bar_OGidnum-gnidnum.png')

OGs.to_csv('out/OG_filter.tsv', sep='\t', index=False)

"""
OUTPUT
Total filtered OGs: 8551

DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_cluster3/connect_OG_graph/connect_OG_graph.py
    ../../ortho_cluster3/connect_OG_graph/out/components.tsv
../../ortho_cluster3/hits2graph/hits2graph.py
    ../../ortho_cluster3/hits2graph/out/hit_graph.tsv
"""