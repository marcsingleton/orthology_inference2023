"""Plot re-aligned alignments."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from src.draw import draw_msa
from src.utils import read_fasta

spid_regex = r'spid=([a-z]+)'

tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}
spids = {tip.name for tip in tree.tips() if tip.name != 'sleb'}

df = pd.read_table('../gap_contrasts/out/total_sums.tsv')
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for i, row in enumerate(head.itertuples()):
        for alignment in ['hmmer', 'mafft']:
            msa = []
            for header, seq in read_fasta(f'../realign_hmmer/out/{alignment}/{row.OGid}.afa'):
                spid = re.search(spid_regex, header).group(1)
                msa.append({'spid': spid, 'seq': seq})
            msa = sorted(msa, key=lambda x: tip_order[x['spid']])  # Re-order sequences

            im = draw_msa([record['seq'].upper() for record in msa])
            plt.imsave(f'out/{label}/{i:03}_{row.OGid}_{alignment}.png', im)

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../gap_contrasts/gap_contrasts.py
    ../gap_contrasts/out/total_sums.tsv
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/hmmer/*.afa
    ../realign_hmmer/out/mafft/*.afa
"""