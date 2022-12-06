"""Draw alignments with largest gap contrasts."""

import os
import re
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.draw import draw_msa
from src.utils import read_fasta

tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}
spids = {tip.name for tip in tree.tips() if tip.name != 'sleb'}
num_contrasts = len(spids) - 1

# 1 PLOT STATISTICS (OGS WITH ALL SPECIES)
df = pd.read_table('out/row_sums.tsv')
df['total'] = df[[f'row{i}' for i in range(num_contrasts)]].sum(axis=1)
contrasts = df.drop(['len1', 'len2', 'total'], axis=1).set_index('OGid').stack()

# 1.1 Tail fraction
max_percentile = 0.3
window_percent = 3
max_idx = ceil(max_percentile * len(df))
window_width = ceil(window_percent / 100 * len(df))

labels = df[[f'row{i}' for i in range(num_contrasts)]].stack().sort_values(ascending=False).index.get_level_values(1)
counts = []
for i in range(num_contrasts):
    count = np.convolve(labels == f'row{i}', np.ones(window_width), 'valid') / window_width
    counts.append((i, count[:max_idx]))
counts = sorted(counts, key=lambda x: x[1].sum(), reverse=True)

x = [i / len(df) for i in range(max_idx)]
for i, count in counts[:9]:
    plt.plot(x, count, label=i, linewidth=1)
plt.plot(x, sum([count for _, count in counts[9:]]), label='others', linewidth=1)
plt.legend()
plt.title(f'Contrast fractions in {window_percent}% sliding windows')
plt.xlabel('Right tail percentile')
plt.ylabel('Fraction')
plt.legend(title='Contrast ID', bbox_to_anchor=(1, 0.5), loc='center left')
plt.subplots_adjust(right=0.8)
plt.savefig('out/line_contrast_window.png')
plt.close()

# 1.2 Overall contrast distribution
plt.hist(contrasts, bins=200)
plt.xlabel('Contrast value')
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
scale = 50
plt.bar(range(num_contrasts), [df[f'row{i}'].mean() for i in range(num_contrasts)],
        yerr=[df[f'row{i}'].std()/scale for i in range(num_contrasts)])
plt.xlabel('Contrast ID')
plt.ylabel(f'Mean ± STD/{scale}')
plt.savefig('out/bar_contrast_mean.png')
plt.close()

# 2 DRAW ALIGNMENTS (ALL OGS)
OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')
df = pd.read_table('out/total_sums.tsv').merge(OG_filter[['OGid', 'ppidnum']], on='OGid', how='left')  # total_sums.tsv has gnidnum already
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.mkdir(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for rank, row in enumerate(head.itertuples()):
        msa = []
        for header, seq in read_fasta(f'../get_repseqs/out/{row.OGid}.afa'):
            spid = re.search(r'spid=([a-z]+)', header).group(1)
            msa.append({'spid': spid, 'seq': seq})
        msa = sorted(msa, key=lambda x: tip_order[x['spid']])

        im = draw_msa([record['seq'] for record in msa])
        plt.imsave(f'out/{label}/{rank:03}_{row.OGid}.png', im)

"""
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../get_repseqs/get_repseqs.py
    ../get_repseqs/out/*.afa
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
./gap_contrasts.py
    ./out/row_sums.tsv
    ./out/total_sums.tsv'
"""