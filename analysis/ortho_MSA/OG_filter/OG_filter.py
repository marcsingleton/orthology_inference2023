"""Filter OGs to representatives by uniqueness criteria and species number."""

import os

import matplotlib.pyplot as plt
import pandas as pd


def get_bitscore(node1, node2, graph):
    """Return bitscore of edge in graph.

    Edges from in-paralogs are not in graph, so they are ignored.
    """
    try:
        return graph[node1][node2]
    except KeyError:
        return 0


def spid_filter(spids):
    conditions = [({'dnov', 'dvir'}, 1),
                  ({'dmoj', 'dnav'}, 1),
                  ({'dinn', 'dgri', 'dhyd'}, 2),
                  ({'dgua', 'dsob'}, 1),
                  ({'dbip', 'dana'}, 1),
                  ({'dser', 'dkik'}, 1),
                  ({'dele', 'drho'}, 1),
                  ({'dtak', 'dbia'}, 1),
                  ({'dsuz', 'dspu'}, 1),
                  ({'dere', 'dtei'}, 1),
                  ({'dsan', 'dyak'}, 1),
                  ({'dmel'}, 1),
                  ({'dmau', 'dsim', 'dsec'}, 1)]
    return all([len(spids & group) >= num for group, num in conditions])


spid_min = 20

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2data[fields['ppid']] = (fields['gnid'], fields['spid'], fields['sqid'])

# Load graph
graph = {}
with open('../../ortho_cluster2/hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        graph[node] = bitscores

# Load GGs
OGid2GGid = {}
with open('../../ortho_cluster2/connect_OG_graph/out/components.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        for OGid in fields['OGids'].split(','):
            OGid2GGid[OGid] = fields['GGid']

# Load OGs
GGid2OGs = {}
with open('../../ortho_cluster2/add_paralogs/out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        component_id, OGid = fields['component_id'], fields['OGid']
        GGid = OGid2GGid[OGid]
        edges = [edge.split(':') for edge in fields['edges'].split(',')]
        ppids = {node for edge in edges for node in edge}
        gnids = {ppid2data[ppid][0] for ppid in ppids}
        spids = {ppid2data[ppid][1] for ppid in ppids}
        sqids = {ppid2data[ppid][2] for ppid in ppids}
        bitscore = 0
        for node1, node2 in edges:
            bitscore += get_bitscore(node1, node2, graph) + get_bitscore(node2, node1, graph)

        if len(gnids) == len(spids) >= spid_min and spid_filter(spids):
            record = {'component_id': component_id, 'OGid': OGid, 'GGid': GGid,
                      'ppidnum': len(ppids), 'sqidnum': len(sqids),
                      'gnidnum': len(gnids), 'spidnum': len(spids),
                      'edgenum': len(edges), 'bitscore': round(bitscore, 1)}
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

counts = OGs['gnidnum'].value_counts()
plt.bar(counts.index, counts.values)
plt.xlabel('Number of single-copy orthologs in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/bar_OGidnum-gnidnum.png')

OGs.to_csv('out/OG_filter.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_cluster2/add_paralogs/add_paralogs.py
    ../../ortho_cluster2/add_paralogs/out/clusters.tsv
../../ortho_cluster2/connect_OG_graph/connect_OG_graph.py
    ../../ortho_cluster2/connect_OG_graph/out/components.tsv
../../ortho_cluster2/hits2graph/hits2graph.py
    ../../ortho_cluster2/hits2graph/out/hit_graph.tsv
"""