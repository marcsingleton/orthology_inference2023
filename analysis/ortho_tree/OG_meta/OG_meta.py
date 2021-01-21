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
with open('../clique4+_community/out/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        spids = set([gnid2spid[gnid] for gnid in gnids])
        sqidnum = sum([gnid2sqidnum[gnid] for gnid in gnids])
        rows.append({'CCid': CCid, 'OGid': OGid,
                     'edgenum': len(edges), 'gnidnum': len(gnids), 'spidnum': len(spids), 'sqidnum': sqidnum})
OGs = pd.DataFrame(rows)

# Print counts
gn27 = OGs['gnidnum'] == 27
sp27 = OGs['spidnum'] == 27
sq27 = OGs['sqidnum'] == 27

print('Total OGs:', len(OGs))
print('OGs with 27 genes:', len(OGs[gn27]))
print('OGs with 27 genes and species:', len(OGs[gn27 & sp27]))
print('OGs with 27 genes, species, and sequences:', len(OGs[gn27 & sp27 & sq27]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_meta.tsv', sep='\t', index=False)

"""
OUTPUT 
Total OGs: 17014
OGs with 27 genes: 7717
OGs with 27 genes and species: 7525
OGs with 27 genes, species, and sequences: 2333

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_community/clique4+_community.py
    ../clique4+_community/out/5clique/gclusters.txt
"""