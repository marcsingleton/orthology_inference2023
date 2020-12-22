"""Extract metadata from OGs including numbers of edges, genes, species, and unique sequences."""

import os
import pandas as pd

# Load gn metadata
gnid2spid = {}
with open('../../ortho_search/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load sqid counts
gnid2sqidnum = {}
with open('../genome_stats/out/gnid_nums.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, gnid, _, sqidnum = line.split()
        gnid2sqidnum[gnid] = int(sqidnum)

# Load OGs
rows = []
with open('../clique5+_community/out/5clique/gclusters.txt') as file:
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

OGs.to_csv('out/OGid2meta.tsv', sep='\t', index=False)

"""
OUTPUT 
Total OGs: 21198
OGs with 27 genes: 7717
OGs with 27 genes and species: 7525
OGs with 27 genes, species, and sequences: 2333

DEPENDENCIES
../../ortho_search/ppid2meta/ppid2meta.py
    ../../ortho_search/ppid2meta/out/ppid2meta.tsv
../clique5+_community/clique5+_community.py
    ../clique5+_community/out/5clique/gclusters.txt
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_nums.tsv
"""