"""Extract metadata from OGs including numbers of edges, genes, species, and sequences."""

import os
import pandas as pd

# Load sequence data
ppid2data, spids = {}, []
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.rstrip('\n').split('\t')
        ppid2data[ppid] = (gnid, spid, sqid)
        spids.append(spid)

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
        GGid, OGids = line.rstrip('\n').split('\t')
        for OGid in OGids.split(','):
            OGid2GGid[OGid] = GGid

# Load OGs
rows = []
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip('\n').split('\t')
        edges = [edge.split(':') for edge in edges.split(',')]
        ppids = {node for edge in edges for node in edge}
        gnids = {ppid2data[ppid][0] for ppid in ppids}
        spids = {ppid2data[ppid][1] for ppid in ppids}
        sqids = {ppid2data[ppid][2] for ppid in ppids}
        bitscore = 0
        for node1, node2 in edges:
            bitscore += graph[node1][node2] + graph[node2][node1]

        row = {'component_id': component_id, 'OGid': OGid, 'GGid': OGid2GGid[OGid],
               'bitscore': round(bitscore, 1), 'edgenum': len(edges),
               'ppidnum': len(ppids), 'sqidnum': len(sqids),
               'gnidnum': len(gnids), 'spidnum': len(spids)}
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
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_cluster3/connect_OG_graph/connect_OG_graph.py
    ../../ortho_cluster3/connect_OG_graph/out/components.tsv
../../ortho_cluster3/hits2graph/hits2graph.py
    ../../ortho_cluster3/hits2graph/out/hit_graph.tsv
"""