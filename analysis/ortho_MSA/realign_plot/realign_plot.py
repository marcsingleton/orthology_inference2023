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


aa2color = {'A': '6dd7a1', 'I': '55c08c', 'L': '55c08c', 'V': '55c08c', 'M': '55c08c',
            'F': 'b897ec', 'Y': 'b897ec', 'W': 'a180d2',
            'S': 'ffbe74', 'T': 'ffbe74',
            'N': '77eaf4', 'Q': '77eaf4',
            'D': 'ee8485', 'E': 'ee8485',
            'H': '96c4ff', 'K': '7fadea', 'R': '7fadea',
            'C': 'faed70', 'G': 'e2dedd', 'P': 'ffb1f1',
            'X': '93908f', '-': 'ffffff', '.': '3f3f3f'}
gap2color = {'-': '3f3f3f', '.': 'ffffff'}

tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
spids = set([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')
df = pd.read_table('../gap_contrasts/out/total_sums.tsv').merge(OG_filter[['OGid', 'sqidnum']], on='OGid', how='left')  # total_sums.tsv has gnidnum already
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

if not os.path.exists('out/norm1/'):
    os.makedirs('out/norm1/')

head = df.sort_values(by='norm1', ascending=False).head(150)
for i, row in enumerate(head.itertuples()):
    msa = load_msa(f'../realign_hmmer/out/{row.OGid}.mfa')

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1].upper() for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa, aa2color=aa2color, gap2color=gap2color)
    plt.imsave(f'out/norm1/{i}_{row.OGid}.png', im)

if not os.path.exists('out/norm2/'):
    os.makedirs('out/norm2/')

head = df.sort_values(by='norm2', ascending=False).head(150)
for i, row in enumerate(head.itertuples()):
    msa = load_msa(f'../realign_hmmer/out/{row.OGid}.mfa')

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1].upper() for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa, aa2color=aa2color, gap2color=gap2color)
    plt.imsave(f'out/norm2/{i}_{row.OGid}.png', im)

"""
DEPENDENCIES
../../../src/draw.py
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
../gap_contrasts/gap_contrasts_calc.py
    ../gap_contrasts/out/total_sums.tsv
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/*.mfa
"""