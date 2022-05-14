"""Extract metadata from OGs including numbers of sequences, genes, species, and edges."""

import os
import pandas as pd


def get_bitscore(node1, node2, graph):
    """Return bitscore of edge in graph.

    Edges from in-paralogs are not in graph, so they are ignored.
    """
    try:
        return graph[node1][node2]
    except KeyError:
        return 0


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
rows = []
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

        row = {'component_id': component_id, 'OGid': OGid, 'GGid': GGid,
               'ppidnum': len(ppids), 'sqidnum': len(sqids),
               'gnidnum': len(gnids), 'spidnum': len(spids),
               'edgenum': len(edges), 'bitscore': round(bitscore, 1)}
        rows.append(row)
OGs = pd.DataFrame(rows)

# Write data to file
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_data.tsv', sep='\t', index=False)

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