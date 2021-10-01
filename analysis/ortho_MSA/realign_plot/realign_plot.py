"""Plot re-aligned alignments."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from src.draw import draw_msa


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                spid = re.search(r'spid=([a-z]+)', line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((spid, seq))
    return msa


tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
spids = set([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')
df = pd.read_table('../gap_contrasts/out/total_sums.tsv').merge(OG_filter[['OGid', 'sqidnum']], on='OGid', how='left')  # total_sums.tsv has gnidnum already
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for i, row in enumerate(head.itertuples()):
        msa = load_msa(f'../realign_hmmer2/out/{row.OGid}.mfa')

        tree = tree_template.shear([spid for spid, _ in msa])
        order = {tip.name: i for i, tip in enumerate(tree.tips())}
        msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
        im = draw_msa(msa)
        plt.imsave(f'out/{label}/{i}_{row.OGid}.png', im)

"""
DEPENDENCIES
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../gap_contrasts/gap_contrasts_calc.py
    ../gap_contrasts/out/total_sums.tsv
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
../realign_hmmer2/realign_hmmer2.py
    ../realign_hmmer2/out/*.mfa
"""