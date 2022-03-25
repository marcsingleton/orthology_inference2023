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

ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'

# Load regions
rows = []
with open('../aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop),
                     'gnidnum': len(ppids.split(',')), 'ppids': ppids.split(',')})
regions = pd.DataFrame(rows)

# Load tree
tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}
spids = {tip.name for tip in tree.tips() if tip.name != 'sleb'}
num_contrasts = len(spids) - 1

# 1 PLOT STATISTICS (OGS WITH ALL SPECIES)
df = pd.read_table('out/row_sums.tsv').merge(regions, on=['OGid', 'start', 'stop'], how='left')
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

# 1.3 Distribution of averages within regions
plt.hist(contrasts.groupby('OGid').mean(), bins=200)
plt.xlabel('Contrast mean in region')
plt.ylabel('Count')
plt.savefig('out/hist_contrast_regionmean.png')
plt.yscale('log')
plt.savefig('out/hist_contrast_regionmean_log.png')
plt.close()

# 1.4 Contrast averages across all regions
plt.bar(list(range(num_contrasts)), [df[f'row{i}'].mean() for i in range(num_contrasts)],
        yerr=[df[f'row{i}'].std()/50 for i in range(num_contrasts)])
plt.xlabel('Contrast ID')
plt.ylabel('Mean ± STD/50')
plt.savefig('out/bar_contrast_mean.png')
plt.close()

# 2 DRAW ALIGNMENTS (ALL OGS)
df = pd.read_table('out/total_sums.tsv').merge(regions[['OGid', 'start', 'stop', 'ppids']], how='left', on=['OGid', 'start', 'stop'])
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.mkdir(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for i, row in enumerate(head.itertuples()):
        msa = []
        for header, seq in read_fasta(f'../insertion_trim/out/{row.OGid}.mfa'):
            ppid = re.search(ppid_regex, header).group(1)
            spid = re.search(spid_regex, header).group(1)
            if ppid in row.ppids:
                msa.append((spid, seq[row.start:row.stop]))

        msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
        im = draw_msa(msa)
        plt.imsave(f'out/{label}/{i:03}_{row.OGid}-{row.start}-{row.stop}.png', im)

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/regions_30.tsv
./contrasts.py
    ./out/row_sums.tsv
    ./out/total_sums.tsv'
"""