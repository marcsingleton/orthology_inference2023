"""Draw test alignments."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from src.draw import draw_msa
from src.utils import read_fasta

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2meta[ppid] = (gnid, sqid)

# Load OGs
rows = []
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        CCid, OGid, _, edges = line.rstrip().split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        for ppid in ppids:
            gnid, sqid = ppid2meta[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'sqid': sqid})
OGs = pd.DataFrame(rows)

# Load pOG metadata and test genes
OG_meta = pd.read_table('../OG_data/out/OG_data.tsv')
genes = pd.read_table('genes.tsv')

# Load tree
tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Draw alignments
if not os.path.exists('out/'):
    os.mkdir('out/')

df = OGs[['gnid', 'OGid']].drop_duplicates().merge(OG_meta, on='OGid', how='right').merge(genes, on='gnid', how='right')
df.to_csv('out/OGs.tsv', sep='\t', index=False)

for row in df.dropna().itertuples():
    if row.sqidnum == row.gnidnum:
        msa = read_fasta(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        msa = read_fasta(f'../align_fastas2/out/{row.OGid}.mfa')
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    msa = [seq for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/{row.OGid}.png', im)

"""
DEPENDENCIES
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.mfa
../OG_data/OG_data.py
    ../OG_data/out/OG_data.tsv
./genes.tsv
"""