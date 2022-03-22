"""Extract metadata from OGs including numbers of edges, genes, species, and sequences."""

import os
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

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
with open('../clique4+_pcommunity/out/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = {node for edge in edges.split('\t') for node in edge.split(',')}
        gnids = {ppid2meta[ppid][0] for ppid in ppids}
        spids = {ppid2meta[ppid][1] for ppid in ppids}
        sqids = {ppid2meta[ppid][2] for ppid in ppids}
        bitscore = 0
        for edge in edges.split('\t'):
            node1, node2 = edge.split(',')
            bitscore += graph[node1][node2] + graph[node2][node1]

        row = {'CCid': CCid, 'OGid': OGid,
               'bitscore': round(bitscore, 1), 'edgenum': len(edges.split('\t')),
               'ppidnum': len(ppids), 'sqidnum': len(sqids),
               'gnidnum': len(gnids), 'spidnum': len(spids)}
        rows.append(row)
OGs = pd.DataFrame(rows)

# Print counts
num = 32
ppnum = OGs['ppidnum'] == num
sqnum = OGs['sqidnum'] == num
gnnum = OGs['gnidnum'] == num
spnum = OGs['spidnum'] == num

print('Total OGs:', len(OGs))
print(f'OGs with {num} genes:', len(OGs[gnnum]))
print(f'OGs with {num} genes and species:', len(OGs[gnnum & spnum]))
print(f'OGs with {num} genes, species, and sequences:', len(OGs[gnnum & spnum & ppnum]))
print(f'OGs with {num} genes, species, and unique sequences:', len(OGs[gnnum & spnum & sqnum]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_meta.tsv', sep='\t', index=False)

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
../clique4+_pcommunity/clique4+_pcommunity.py
    ../../ortho_cluster3/clique4+_pcommunity/out/4clique/pclusters.txt
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""