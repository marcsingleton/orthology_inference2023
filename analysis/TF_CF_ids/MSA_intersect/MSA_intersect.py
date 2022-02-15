"""Intersect TF and CF lists with curated MSAs."""

import os
import pandas as pd

# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2gnid[ppid] = gnid

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            gnid = ppid2gnid[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid})
OG_gnids = pd.DataFrame(rows).drop_duplicates()  # All OGs with genes
OG_filter = pd.read_table('../../ortho_MSA/OG_filter/out/OG_filter.tsv', usecols=['CCid', 'OGid', 'gOGid'])  # OGs after filtering
OGs = OG_filter.merge(OG_gnids, how='left', on=['OGid', 'CCid'])  # Filtered OGs with genes

# Load TFs and CFs and merge with OGs
TFs = pd.read_table('../update_ids/out/TFs.txt', names=['gnid'])
CFs = pd.read_table('../update_ids/out/CFs.txt', names=['gnid'])

OGs_TF = TFs.merge(OGs[['OGid', 'gnid']], how='inner', on=['gnid'])
OGs_CF = CFs.merge(OGs[['OGid', 'gnid']], how='inner', on=['gnid'])

# Print stats and write to file
print('Number of TF OGs:', len(OGs_TF))
print('Number of unique TF genes:', OGs_TF['gnid'].nunique())
print('Number of CF OGs:', len(OGs_CF))
print('Number of unique CF genes:', OGs_CF['gnid'].nunique())

if not os.path.exists('out/'):
    os.mkdir('out/')

OGs_TF.to_csv('out/TFs.tsv', sep='\t', index=False)
OGs_CF.to_csv('out/CFs.tsv', sep='\t', index=False)

"""
OUTPUT
Number of TF OGs: 563
Number of unique TF genes: 563
Number of CF OGs: 251
Number of unique CF genes: 251

NOTES
The number of unique genes for each set of OGs is the same as the total. This means the genes are uniquely mapped to
OGs (as opposed to one gene potentially being in multiple OGs which is quite possible). The filtered list of OGs already
removed OGs with multiple representatives from each species, so the OGs also uniquely map to genes.

DEPENDENCIES
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_MSA/OG_filter/OG_filter.py
    ../../ortho_MSA/OG_filter/out/OG_filter.tsv
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../update_ids/update_ids.py
    ../update_ids/out/*.txt
"""