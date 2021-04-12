"""Draw alignments with largest gap contrasts."""

import os
import re
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
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


OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')
tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
spids = set([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])
num_contrasts = len(spids) - 1

# 1 PLOT STATISTICS (OGS WITH ALL SPECIES)
df = pd.read_table('out/row_sums.tsv').merge(OG_filter, on='OGid', how='left')
df['total'] = df[[f'row{i}' for i in range(num_contrasts)]].sum(axis=1)
df['avg'] = df['total'] / df['len2']

# 1.1 Tail fraction
max_percentile = 0.30
window_percent = 3
window_width = ceil(window_percent / 200 * len(df))
max_idx = int(max_percentile * len(df)) - 1  # Subtract 1 due to 0 indexing
x = [i / len(df) for i in range(window_width + 1, max_idx - window_width + 1)]

contrasts = pd.concat([df[f'row{i}'] for i in range(num_contrasts)], keys=[f'row{i}' for i in range(num_contrasts)],
                      names=['contrast_id', 'OGid']).sort_values(ascending=False)
counts = [(i, contrasts[:max_idx].index.get_level_values('contrast_id') == f'row{i}') for i in range(num_contrasts)]
tails = sorted([(i, np.convolve(count, np.ones(2 * window_width + 1) / (2 * window_width + 1), 'valid')) for i, count in counts],
               key=lambda y: sum(y[1]), reverse=True)
for i, tail in tails[:9]:
    plt.plot(x, tail, label=i, linewidth=1)
plt.plot(x, sum([tail for _, tail in tails[9:]]), label='others', linewidth=1)
plt.legend()
plt.title(f'Contrast Fractions in {window_percent}% Sliding Windows')
plt.xlabel('Right Tail Percentile')
plt.ylabel('Fraction')
plt.legend(title='Contrast ID', bbox_to_anchor=(1.025, 0.5), loc='center left')
plt.subplots_adjust(right=0.8)
plt.savefig('out/line_contrast_window.png')
plt.close()

# 1.2 Overall contrast distribution
plt.hist(contrasts, bins=200)
plt.xlabel('Contrast Value')
plt.ylabel('Count')
plt.savefig('out/hist_contrast.png')
plt.yscale('log')
plt.savefig('out/hist_contrast_log.png')
plt.close()

# 1.3 Distribution of averages within OGs
plt.hist(contrasts.groupby('OGid').mean(), bins=200)
plt.xlabel('Contrast mean in OG')
plt.ylabel('Count')
plt.savefig('out/hist_contrast_OGmean.png')
plt.yscale('log')
plt.savefig('out/hist_contrast_OGmean_log.png')
plt.close()

# 1.4 Contrast averages across all OGs
plt.bar(list(range(num_contrasts)), [df[f'row{i}'].mean() for i in range(num_contrasts)],
        yerr=[df[f'row{i}'].std()/50 for i in range(num_contrasts)])
plt.xlabel('Contrast ID')
plt.ylabel('Mean ± STD/50')
plt.savefig('out/bar_contrast_mean.png')
plt.close()

# 2 DRAW ALIGNMENTS (ALL OGS)
df = pd.read_table('out/total_sums.tsv').merge(OG_filter[['OGid', 'sqidnum']], on='OGid', how='left')  # total_sums.tsv has gnidnum already
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

if not os.path.exists('out/norm1/'):
    os.mkdir('out/norm1/')

# 2.1 Ranked by sum
head1 = df.sort_values(by='norm1', ascending=False).head(150)
for i, row in enumerate(head1.itertuples()):
    if row.sqidnum == row.gnidnum:
        msa = load_msa(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        msa = load_msa(f'../align_fastas2-2/out/{row.OGid}.mfa')

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1] for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/norm1/{i}_{row.OGid}.png', im)

# 2.2 Ranked by avg
if not os.path.exists('out/norm2/'):
    os.mkdir('out/norm2/')

head1 = df.sort_values(by='norm2', ascending=False).head(150)
for i, row in enumerate((head1.itertuples())):
    if row.sqidnum == row.gnidnum:
        msa = load_msa(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        msa = load_msa(f'../align_fastas2-2/out/{row.OGid}.mfa')

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1] for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/norm2/{i}_{row.OGid}.png', im)

"""
../../../src/draw.py
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
./gap_contrasts_calc.py
    ./out/row_sums.tsv
"""