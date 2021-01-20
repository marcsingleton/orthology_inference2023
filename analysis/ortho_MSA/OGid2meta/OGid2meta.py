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
with open('../../ortho_cluster3/genome_stats/out/gnid_nums.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, gnid, _, sqidnum = line.split()
        gnid2sqidnum[gnid] = int(sqidnum)

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        spids = set([gnid2spid[gnid] for gnid in gnids])
        sqidnum = sum([gnid2sqidnum[gnid] for gnid in gnids])
        rows.append({'CCid': CCid, 'OGid': OGid,
                     'edgenum': len(edges), 'gnidnum': len(gnids), 'spidnum': len(spids), 'sqidnum': sqidnum})
OGs = pd.DataFrame(rows)

# Print counts
gn26 = OGs['gnidnum'] == 26
sp26 = OGs['spidnum'] == 26
sq26 = OGs['sqidnum'] == 26

print('Total OGs:', len(OGs))
print('OGs with 26 genes:', len(OGs[gn26]))
print('OGs with 26 genes and species:', len(OGs[gn26 & sp26]))
print('OGs with 26 genes, species, and sequences:', len(OGs[gn26 & sp26 & sq26]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OGid2meta.tsv', sep='\t', index=False)

"""
OUTPUT 
Total OGs: 16757
OGs with 26 genes: 7904
OGs with 26 genes and species: 7733
OGs with 26 genes, species, and sequences: 2492

DEPENDENCIES
../../ortho_cluster3/genome_stats/genome_stats.py
    ../../ortho_cluster3/genome_stats/out/gnid_nums.tsv
../../ortho_search/ppid2meta/ppid2meta.py
    ../../ortho_search/ppid2meta/out/ppid2meta.tsv
../../ortho_cluster3/clique4+_community/clique4+_community.py
    ../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt
"""