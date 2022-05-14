"""Draw test alignments."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from src.draw import draw_msa
from src.utils import read_fasta

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2gnid[fields['ppid']] = fields['gnid']

# Load OGs
rows = []
with open('../../ortho_cluster2/add_paralogs/out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
        for ppid in ppids:
            rows.append({'component_id': fields['component_id'], 'OGid': fields['OGid'], 'ppid': ppid, 'gnid': ppid2gnid[ppid]})
OGs = pd.DataFrame(rows)

# Load OG data and test genes
OG_data = pd.read_table('../OG_data/out/OG_data.tsv')
genes = pd.read_table('genes.tsv')

# Load tree
tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Draw alignments
if not os.path.exists('out/'):
    os.mkdir('out/')

df = OGs[['gnid', 'OGid']].drop_duplicates().merge(OG_data, on='OGid', how='right').merge(genes, on='gnid', how='right')
df.to_csv('out/OGs.tsv', sep='\t', index=False)

for row in df.dropna().itertuples():
    if row.ppidnum == row.gnidnum:
        msa = read_fasta(f'../align_fastas1/out/{row.OGid}.afa')
    else:
        msa = read_fasta(f'../align_fastas2/out/{row.OGid}.afa')
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    msa = [seq for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/{row.OGid}.png', im)

"""
DEPENDENCIES
../../ortho_cluster2/add_paralogs/add_paralogs.py
    ../../ortho_cluster2/add_paralogs/out/clusters.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.afa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.afa
../OG_data/OG_data.py
    ../OG_data/out/OG_data.tsv
./genes.tsv
"""