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
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2meta[ppid] = (gnid, sqid)

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            gnid, sqid = ppid2meta[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'sqid': sqid})
OGs = pd.DataFrame(rows)

# Load pOG metadata and test genes
OG_meta = pd.read_table('../OG_meta/out/OG_meta.tsv')
genes = pd.read_table('genes.tsv')

# Load tree
tree = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
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
        msa = read_fasta(f'../align_fastas2-2/out/{row.OGid}.mfa')
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    msa = [seq for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/{row.OGid}.png', im)

"""
DEPENDENCIES
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
./genes.tsv
"""