"""Draw alignments with largest gap contrasts."""

import os
import re
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.draw import draw_alignment


def load_alignment(path):
    MSA = []
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
            MSA.append((spid, seq))
    return MSA


tree = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree = tree.shear([tip.name for tip in tree.tips() if tip.name != 'sleb'])

OGid2meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv').drop(['CCid', 'edgenum'], axis=1)
df = pd.read_table('out/row_sums.tsv').merge(OGid2meta, on='OGid', how='left')
df['c_sum'] = df[[f'c{i}' for i in range(25)]].sum(axis=1)
df['c_avg'] = df['c_sum'] / df['len2']

# 1 PLOT STATISTICS
# 1.1 Tail fraction
num_contrasts = 25
max_percentile = 0.30

size = ceil(0.015 * len(df))
max_idx = int(max_percentile * len(df)) - 1  # Subtract 1 due to 0 indexing
x = [i / len(df) for i in range(size + 1, max_idx - size + 1)]

contrasts = pd.concat([df[f'c{i}'] for i in range(num_contrasts)], keys=[f'c{i}' for i in range(num_contrasts)],
                      names=['contrast_id', 'OGid']).sort_values(ascending=False)
counts = [(i, contrasts[:max_idx].index.get_level_values('contrast_id') == f'c{i}') for i in range(num_contrasts)]
tails = sorted([(i, np.convolve(count, np.ones(2 * size + 1) / (2 * size + 1), 'valid')) for i, count in counts],
               key=lambda y: sum(y[1]), reverse=True)
for i, tail in tails[:9]:
    plt.plot(x, tail, label=f'c{i}', linewidth=1)
plt.plot(x, sum([tail for _, tail in tails[9:]]), label='others', linewidth=1)
plt.legend()
plt.title('Contrast Fraction in 3% Sliding Windows')
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
plt.bar(list(range(num_contrasts)), [df[f'c{i}'].mean() for i in range(num_contrasts)],
        yerr=[df[f'c{i}'].std()/50 for i in range(num_contrasts)])
plt.xlabel('Contrast ID')
plt.ylabel('Mean ± STD/50')
plt.savefig('out/bar_contrast_mean.png')
plt.close()

# 2 DRAW ALIGNMENTS
if not os.path.exists('out/sum/'):
    os.mkdir('out/sum/')

# 2.1 Ranked by sum
head1 = df.sort_values(by='c_sum', ascending=False).head(100)
for i, row in enumerate((head1.itertuples())):
    if row.sqidnum == 26:
        MSA = load_alignment(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-1/out/{row.OGid}.mfa')

    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(MSA, f'out/sum/{i}_{row.OGid}.png')

# 2.2 Ranked by avg
if not os.path.exists('out/avg/'):
    os.mkdir('out/avg/')

head1 = df.sort_values(by='c_avg', ascending=False).head(100)
for i, row in enumerate((head1.itertuples())):
    if row.sqidnum == 26:
        MSA = load_alignment(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-2/out/{row.OGid}.mfa')

    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(MSA, f'out/avg/{i}_{row.OGid}.png')

"""
../../../src/draw.py
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
./align_contrasts_calc.py
    ./out/row_sums.tsv
"""