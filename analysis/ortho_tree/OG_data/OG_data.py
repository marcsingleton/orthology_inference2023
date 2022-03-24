"""Extract metadata from OGs including numbers of edges, genes, species, and sequences."""

import os
import pandas as pd

# Load seq metadata
ppid2meta, spids = {}, []
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)
        spids.append(spid)

# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        graph[node] = bitscores

# Load OGs
rows = []
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip().split('\t')
        edges = [edge.split(':') for edge in edges.split(',')]
        ppids = {node for edge in edges for node in edge}
        gnids = {ppid2meta[ppid][0] for ppid in ppids}
        spids = {ppid2meta[ppid][1] for ppid in ppids}
        sqids = {ppid2meta[ppid][2] for ppid in ppids}
        bitscore = 0
        for node1, node2 in edges:
            bitscore += graph[node1][node2] + graph[node2][node1]

        row = {'component_id': component_id, 'OGid': OGid,
               'bitscore': round(bitscore, 1), 'edgenum': len(edges),
               'ppidnum': len(ppids), 'sqidnum': len(sqids),
               'gnidnum': len(gnids), 'spidnum': len(spids)}
        rows.append(row)
OGs = pd.DataFrame(rows)

# Print counts
spidnum = len(spids)
ppid_filter = OGs['ppidnum'] == spidnum
gnid_filter = OGs['gnidnum'] == spidnum
spid_filter = OGs['spidnum'] == spidnum

print('Total OGs:', len(OGs))
print(f'OGs with {spidnum} genes:', len(OGs[gnid_filter]))
print(f'OGs with {spidnum} genes and species:', len(OGs[gnid_filter & spid_filter]))
print(f'OGs with {spidnum} genes, species, and sequences:', len(OGs[gnid_filter & spid_filter & ppid_filter]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_data.tsv', sep='\t', index=False)

"""
OUTPUT 
Total OGs: 23242
OGs with 32 genes: 7495
OGs with 32 genes and species: 7267
OGs with 32 genes, species, and sequences: 1919
OGs with 32 genes, species, and unique sequences: 4864

DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""