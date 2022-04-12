"""Intersect TF and CF lists with curated MSAs."""

import os
import pandas as pd

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.rstrip('\n').split('\t')
        ppid2gnid[ppid] = gnid

# Load OGs
rows = []
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip('\n').split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        for ppid in ppids:
            gnid = ppid2gnid[ppid]
            rows.append({'component_id': component_id, 'OGid': OGid, 'gnid': gnid})
OG_gnids = pd.DataFrame(rows).drop_duplicates()  # All OGs with genes
OG_filter = pd.read_table('../../ortho_MSA/OG_filter/out/OG_filter.tsv', usecols=['component_id', 'OGid', 'GGid'])  # OGs after filtering
OGs = OG_filter.merge(OG_gnids, how='left', on=['OGid', 'component_id'])  # Filtered OGs with genes

# Load TFs and CFs and merge with OGs
TFs = pd.read_table('../update_ids/out/TFs.txt', names=['gnid'])
CFs = pd.read_table('../update_ids/out/CFs.txt', names=['gnid'])

OGs_TF = TFs.merge(OGs[['OGid', 'gnid']], how='inner', on=['gnid'])
OGs_CF = CFs.merge(OGs[['OGid', 'gnid']], how='inner', on=['gnid'])

# Write stats and IDs to file
if not os.path.exists('out/'):
    os.mkdir('out/')

output = f"""\
Number of TF OGs: {len(OGs_TF)}
Number of unique TF genes: {OGs_TF['gnid'].nunique()}
Number of CF OGs: {len(OGs_CF)}
Number of unique CF genes: {OGs_CF['gnid'].nunique()}
"""
with open('out/output.txt', 'w') as file:
    file.write(output)

OGs_TF.to_csv('out/TFs.tsv', sep='\t', index=False)
OGs_CF.to_csv('out/CFs.tsv', sep='\t', index=False)

"""
NOTES
The number of unique genes for each set of OGs is the same as the total. This means the genes are uniquely mapped to
OGs (as opposed to one gene potentially being in multiple OGs which is quite possible). The filtered list of OGs already
removed OGs with multiple representatives from each species, so the OGs also uniquely map to genes.

DEPENDENCIES
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_MSA/OG_filter/OG_filter.py
    ../../ortho_MSA/OG_filter/out/OG_filter.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../update_ids/update_ids.py
    ../update_ids/out/*.txt
"""