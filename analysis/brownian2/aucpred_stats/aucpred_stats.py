"""Plot statistics of filtered OGs."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace
from sklearn.decomposition import PCA


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


ppid_regex = r'ppid=([A-Za-z0-9_]+)'
min_lengths = sorted([int(path.split('_')[1][:-4]) for path in os.listdir('../aucpred_filter/out/') if path.endswith('.tsv')])
pdidx = pd.IndexSlice

# Load features
features = pd.read_table('../get_features/out/features.tsv')
features.loc[features['kappa'] == -1, 'kappa'] = 1
features.loc[features['omega'] == -1, 'omega'] = 1

# Parse regions
rows = []
for min_length in min_lengths:
    with open(f'../aucpred_filter/out/regions_{min_length}.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            OGid, start, stop, disorder, ppids = line.split()

            msa = load_msa(f'../insertion_trim/out/{OGid}.mfa')
            msa = {re.search(ppid_regex, header).group(1): seq for header, seq in msa}

            for ppid in ppids.split(','):
                segment = msa[ppid][int(start):int(stop)]
                length = len([sym for sym in segment if sym not in ['-', '.']])
                rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'disorder': disorder == 'True',
                             'ppid': ppid, 'min_length': min_length, 'length': length})
df = pd.DataFrame(rows)

# Plots of combined segment sets
if not os.path.exists('out/'):
    os.mkdir('out/')

# Number of regions by length threshold
disorder, order = [], []
for min_length in min_lengths:
    disorder.append(len(df.loc[df['disorder'] & (df['min_length'] == min_length), ['OGid', 'start', 'stop']].drop_duplicates()))
    order.append(len(df.loc[~df['disorder'] & (df['min_length'] == min_length), ['OGid', 'start', 'stop']].drop_duplicates()))
plt.plot(min_lengths, disorder, color='C0', label='disorder')
plt.plot(min_lengths, order, color='C1', label='order')
plt.xlabel('Length threshold')
plt.ylabel('Number of regions')
plt.legend()
plt.savefig('out/line_numregions-minlength.png')
plt.close()

# Number of OGs by length threshold
disorder, order = [], []
for min_length in min_lengths:
    disorder.append(len(df.loc[df['disorder'] & (df['min_length'] == min_length), 'OGid'].drop_duplicates()))
    order.append(len(df.loc[~df['disorder'] & (df['min_length'] == min_length), 'OGid'].drop_duplicates()))
plt.plot(min_lengths, disorder, color='C0', label='disorder')
plt.plot(min_lengths, order, color='C1', label='order')
plt.xlabel('Length threshold')
plt.ylabel('Number of unique OGs')
plt.legend()
plt.savefig('out/line_numOGs-minlength.png')
plt.close()

# Plots of individual segment sets
for min_length in min_lengths:
    segments = df[df['min_length'] == min_length].merge(features, how='left', on=['OGid', 'start', 'stop', 'ppid'])
    regions = segments.groupby(['OGid', 'start', 'stop', 'disorder'])
    mean = regions.mean()

    if not os.path.exists(f'out/regions_{min_length}/'):
        os.mkdir(f'out/regions_{min_length}/')

    # Mean region length histogram
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = mean['length'].min(), mean['length'].max()
    axs[0].hist(mean.loc[pdidx[:, :, :, True], 'length'], bins=linspace(xmin, xmax, 100), color='C0', label='disorder')
    axs[1].hist(mean.loc[pdidx[:, :, :, False], 'length'], bins=linspace(xmin, xmax, 100), color='C1', label='order')
    axs[1].set_xlabel('Mean length of region')
    for i in range(2):
        axs[i].set_ylabel('Number of regions')
        axs[i].legend()
    plt.savefig(f'out/regions_{min_length}/hist_numregions-length.png')
    plt.close()

    # Number of sequences in region bar plot
    fig, ax = plt.subplots()
    counts1 = regions.size()[pdidx[:, :, :, True]].value_counts()
    counts2 = regions.size()[pdidx[:, :, :, False]].value_counts()
    ax.bar(counts1.index - 0.35/2, counts1.values, color='C0', label='disorder', width=0.35)
    ax.bar(counts2.index + 0.35/2, counts2.values, color='C1', label='order', width=0.35)
    ax.set_xlabel('Number of sequences in region')
    ax.set_ylabel('Number of regions')
    ax.legend()
    plt.savefig(f'out/regions_{min_length}/bar_numregions-numseqs.png')
    plt.close()

    # Counts of regions and unique OGs in each class
    disorder = segments[segments['disorder']]
    order = segments[~segments['disorder']]

    plt.bar([0, 1], [len(disorder[['OGid', 'start', 'stop']].drop_duplicates()), len(order[['OGid', 'start', 'stop']].drop_duplicates())],
            tick_label=['disorder', 'order'], color=['C0', 'C1'], width=0.35)
    plt.xlim((-0.5, 1.5))
    plt.ylabel('Number of regions')
    plt.savefig(f'out/regions_{min_length}/bar_numregions-DO.png')
    plt.close()

    plt.bar([0, 1], [len(disorder['OGid'].drop_duplicates()), len(order['OGid'].drop_duplicates())],
            tick_label=['disorder', 'order'], color=['C0', 'C1'], width=0.35)
    plt.xlim((-0.5, 1.5))
    plt.ylabel('Number of unique OGs')
    plt.savefig(f'out/regions_{min_length}/bar_numOGs-DO.png')
    plt.close()

    # Feature variance pie chart
    var = mean.var().drop(['min_length', 'length']).sort_values(ascending=False)  # Remove "non-feature" columns
    truncate = pd.concat([var[:4], pd.Series({'other': var[4:].sum()})])
    plt.pie(truncate.values, labels=truncate.index)
    plt.title(f'Feature variance, length ≥ {min_length}')
    plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(right=0.7)
    plt.savefig(f'out/regions_{min_length}/pie_variance.png')
    plt.close()

    # Feature PCAs
    pca = PCA(n_components=5)
    idx = mean.index.get_level_values('disorder').array.astype(bool)
    x = mean.drop(['min_length', 'length'], axis=1)  # Remove "non-feature" columns

    transform = pca.fit_transform(x.to_numpy())
    plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, alpha=0.05, edgecolors='none')
    plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, alpha=0.05, edgecolors='none')
    plt.title(f'unnormed, length ≥ {min_length}')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    legend = plt.legend(markerscale=2)
    for lh in legend.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/regions_{min_length}/pca_unnorm.png')
    plt.close()

    norm = (x - x.mean()) / x.std()
    transform = pca.fit_transform(norm.to_numpy())
    plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, alpha=0.05, edgecolors='none')
    plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, alpha=0.05, edgecolors='none')
    plt.title(f'z-score, length ≥ {min_length}')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    legend = plt.legend(markerscale=2)
    for lh in legend.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/regions_{min_length}/pca_z-score.png')
    plt.close()

    norm = (x - x.min()) / (x.max()-x.min())
    transform = pca.fit_transform(norm.to_numpy())
    plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, alpha=0.05, edgecolors='none')
    plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, alpha=0.05, edgecolors='none')
    plt.title(f'min-max, length ≥ {min_length}')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    legend = plt.legend(markerscale=2)
    for lh in legend.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/regions_{min_length}/pca_min-max.png')
    plt.close()

"""
DEPENDENCIES
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/segments_*.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""