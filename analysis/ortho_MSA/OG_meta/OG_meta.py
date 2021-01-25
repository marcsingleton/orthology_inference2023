"""Extract metadata from OGs including numbers of edges, genes, species, and unique sequences."""

import os
import pandas as pd

# Load seq metadata
gnid2spid = {}
gnid2sqidnum = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        _, gnid, spid, repr = line.split()
        gnid2spid[gnid] = spid
        if repr == 'True':
            gnid2sqidnum[gnid] = gnid2sqidnum.get(gnid, 0) + 1

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        spids = set([gnid2spid[gnid] for gnid in gnids])
        sqidnum = sum([gnid2sqidnum[gnid] for gnid in gnids])
        rows.append({'CCid': CCid, 'OGid': OGid,
                     'edgenum': len(edges.split('\t')), 'gnidnum': len(gnids), 'spidnum': len(spids), 'sqidnum': sqidnum})
OGs = pd.DataFrame(rows)

# Print counts
num = 26
gnnum = OGs['gnidnum'] == num
spnum = OGs['spidnum'] == num
sqnum = OGs['sqidnum'] == num

print('Total OGs:', len(OGs))
print(f'OGs with {num} genes:', len(OGs[gnnum]))
print(f'OGs with {num} genes and species:', len(OGs[gnnum & spnum]))
print(f'OGs with {num} genes, species, and sequences:', len(OGs[gnnum & spnum & sqnum]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_meta.tsv', sep='\t', index=False)

"""
OUTPUT 
Total OGs: 14533
OGs with 26 genes: 7903
OGs with 26 genes and species: 7733
OGs with 26 genes, species, and sequences: 2492

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/clique4+_gcommunity/clique4+_gcommunity2.py
    ../../ortho_cluster3/clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt
"""