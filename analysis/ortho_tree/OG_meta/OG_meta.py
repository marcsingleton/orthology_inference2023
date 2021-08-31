"""Extract metadata from OGs including numbers of edges, genes, species, and sequences."""

import os
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

# Load pgraph
pgraph = {}
with open('../hits2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        pgraph[node] = bitscores

# Load OGs
rows = []
with open('../clique4+_pcommunity/out/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        gnids = set([ppid2meta[ppid][0] for ppid in ppids])
        spids = set([ppid2meta[ppid][1] for ppid in ppids])
        sqids = set([ppid2meta[ppid][2] for ppid in ppids])
        bitscore = 0
        for edge in edges.split('\t'):
            node1, node2 = edge.split(',')
            bitscore += pgraph[node1][node2] + pgraph[node2][node1]

        d = {'CCid': CCid, 'OGid': OGid,
             'bitscore': round(bitscore, 1), 'edgenum': len(edges.split('\t')),
             'ppidnum': len(ppids), 'sqidnum': len(sqids),
             'gnidnum': len(gnids), 'spidnum': len(spids)}
        rows.append(d)
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
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_pcommunity/clique4+_pcommunity.py
    ../../ortho_cluster3/clique4+_pcommunity/out/4clique/pclusters.txt
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph.tsv
"""