"""Plot statistics associated with contrasts."""

import os
import re
from math import exp, pi

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from numpy import linspace
from src.draw import plot_msa
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


spid_regex = r'spid=([a-z]+)'
pdidx = pd.IndexSlice

# Load contrasts and tree
contrasts = pd.read_table('../get_contrasts/out/contrasts.tsv')
tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)

# Load seq metadata
ppid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, _, spid, _ = line.split()
        ppid2spid[ppid] = spid

# Parse segments
rows, region2spids = [], {}
with open('../aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'disorder': disorder == 'True'})
        region2spids[(OGid, int(start), int(stop))] = [ppid2spid[ppid] for ppid in ppids.split(',')]
segments = pd.DataFrame(rows)
df1 = segments.merge(contrasts, how='right', on=['OGid', 'start', 'stop']).set_index(['OGid', 'start', 'stop', 'disorder', 'contrast_id'])

# 1 CONTRASTS
if not os.path.exists('out/contrasts/'):
    os.makedirs('out/contrasts/')

# 1.1 Plot contrast distributions
disorder = df1.loc[pdidx[:, :, :, True, :], :]
order = df1.loc[pdidx[:, :, :, False, :], :]
for feature_label in df1.columns:
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = df1[feature_label].min(), df1[feature_label].max()
    axs[0].hist(disorder[feature_label], bins=linspace(xmin, xmax, 150), color='C0', label='disorder')
    axs[1].hist(order[feature_label], bins=linspace(xmin, xmax, 150), color='C1', label='order')
    axs[1].set_xlabel(f'Contrast value ({feature_label})')
    for i in range(2):
        axs[i].set_ylabel('Number of contrasts')
        axs[i].legend()
    plt.savefig(f'out/contrasts/{feature_label}.png')
    for i in range(2):
        axs[i].set_yscale('log')
    plt.savefig(f'out/contrasts/{feature_label}_log.png')
    plt.close()

# 1.2 Plot regions with extreme contrasts
df2 = df1.abs()
for feature_label in df2.columns:
    if not os.path.exists(f'out/contrasts/{feature_label}/'):
        os.makedirs(f'out/contrasts/{feature_label}/')

    regions = set()
    ranked = df2[feature_label].sort_values(ascending=False).iteritems()
    while len(regions) < 20:
        (OGid, start, stop, _, _), _ = next(ranked)
        if (OGid, start, stop) in regions:
            continue
        regions.add((OGid, start, stop))

        # Load MSA, filter seqs, and re-order
        msa1 = load_msa(f'../insertion_trim/out/{OGid}.mfa')
        msa1 = {re.search(spid_regex, header).group(1): seq for header, seq in msa1}

        spids = region2spids[(OGid, start, stop)]
        tree = tree_template.shear(spids)
        order = {tip.name: i for i, tip in enumerate(tree.tips())}
        msa2 = [msa1[spid].upper()[start:stop] for spid in sorted(spids, key=lambda x: order[x])]
        fig = plot_msa(msa2, figsize=(8, 6), x_start=start)
        plt.savefig(f'out/contrasts/{feature_label}/{len(regions)-1}_{OGid}-{start}-{stop}.png', bbox_inches='tight')
        plt.close()

# 2 CONTRAST MEANS
if not os.path.exists('out/means/'):
    os.makedirs('out/means/')

# 2.1 Plot contrast means distributions
mean = df1.groupby(['OGid', 'start', 'stop', 'disorder']).mean()
disorder = mean.loc[pdidx[:, :, :, True, :], :]
order = mean.loc[pdidx[:, :, :, False, :], :]
for feature_label in mean.columns:
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = mean[feature_label].min(), mean[feature_label].max()
    axs[0].hist(disorder[feature_label], bins=linspace(xmin, xmax, 150), color='C0', label='disorder')
    axs[1].hist(order[feature_label], bins=linspace(xmin, xmax, 150), color='C1', label='order')
    axs[1].set_xlabel(f'Mean contrast value ({feature_label})')
    for i in range(2):
        axs[i].set_ylabel('Number of regions')
        axs[i].legend()
    plt.savefig(f'out/means/{feature_label}.png')
    plt.close()

# 2.2 Plot standardized contrast means distributions
# These are sample means which have a mean of 0 and variance sigma^2/n
# Estimate sigma^2 from contrasts by mean of contrast squares (since theoretical contrast mean is 0)
# Standardize sample means by dividing by sigma/sqrt(n)
# Regions with constant contrasts will have 0 variance, so the normalization will result in a NaN
# While it is possible for a tree with unequal tip values to yield constant (non-zero) contrasts, it is unlikely
# Thus constant contrasts are assumed to equal zero
var = ((df1**2).groupby(['OGid', 'start', 'stop', 'disorder']).mean())
size = df1.groupby(['OGid', 'start', 'stop', 'disorder']).size()
std = (mean / (var.div(size, axis=0))**0.5).fillna(0)
disorder = std.loc[pdidx[:, :, :, True, :], :]
order = std.loc[pdidx[:, :, :, False, :], :]
for feature_label in std.columns:
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = std[feature_label].min(), std[feature_label].max()
    axs[0].hist(disorder[feature_label], density=True, bins=linspace(xmin, xmax, 150), color='C0', label='disorder')
    axs[1].hist(order[feature_label], density=True, bins=linspace(xmin, xmax, 150), color='C1', label='order')
    axs[1].set_xlabel(f'Standardized mean contrast value ({feature_label})')
    for i in range(2):
        axs[i].plot(linspace(xmin, xmax), [1/(2*pi)**0.5 * exp(-x**2/2) for x in linspace(xmin, xmax)], color='black')
        axs[i].set_ylabel('Density of regions')
        axs[i].legend()
    plt.savefig(f'out/means/{feature_label}_std.png')
    plt.close()

# 2.3 Plot contrast mean PCAs
pca = PCA(n_components=10)
idx = mean.index.get_level_values('disorder').array.astype(bool)

transform = pca.fit_transform(mean.to_numpy())
plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, color='C0', alpha=0.1, edgecolors='none')
plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, color='C1', alpha=0.1, edgecolors='none')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('raw means')
legend = plt.legend(markerscale=2)
for lh in legend.legendHandles:
    lh.set_alpha(1)
plt.savefig(f'out/means/pca_mean.png')
plt.close()

idx = std.index.get_level_values('disorder').array.astype(bool)
transform = pca.fit_transform(std.to_numpy())
plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, color='C0', alpha=0.1, edgecolors='none')
plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, color='C1', alpha=0.1, edgecolors='none')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('standardized means')
legend = plt.legend(markerscale=2)
for lh in legend.legendHandles:
    lh.set_alpha(1)
plt.savefig(f'out/means/pca_mean_std.png')
plt.close()

# 2.4 Plot regions with extreme means
df2 = mean.abs()
for feature_label in df2.columns:
    if not os.path.exists(f'out/means/{feature_label}/'):
        os.makedirs(f'out/means/{feature_label}/')

    regions = set()
    ranked = df2[feature_label].sort_values(ascending=False).iteritems()
    while len(regions) < 20:
        (OGid, start, stop, _), _ = next(ranked)
        if (OGid, start, stop) in regions:
            continue
        regions.add((OGid, start, stop))

        # Load MSA, filter seqs, and re-order
        msa1 = load_msa(f'../insertion_trim/out/{OGid}.mfa')
        msa1 = {re.search(spid_regex, header).group(1): seq for header, seq in msa1}

        spids = region2spids[(OGid, start, stop)]
        tree = tree_template.shear(spids)
        order = {tip.name: i for i, tip in enumerate(tree.tips())}
        msa2 = [msa1[spid].upper()[start:stop] for spid in sorted(spids, key=lambda x: order[x])]
        fig = plot_msa(msa2, figsize=(8, 6), x_start=start)
        plt.savefig(f'out/means/{feature_label}/{len(regions)-1}_{OGid}-{start}-{stop}.png', bbox_inches='tight')
        plt.close()

# 3 RATES
if not os.path.exists('out/rates/'):
    os.makedirs('out/rates/')

# 3.1 Plot rate distributions
rate = ((df1**2).groupby(['OGid', 'start', 'stop', 'disorder']).mean())
disorder = rate.loc[pdidx[:, :, :, True, :], :]
order = rate.loc[pdidx[:, :, :, False, :], :]
for feature_label in rate.columns:
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = rate[feature_label].min(), rate[feature_label].max()
    axs[0].hist(disorder[feature_label], bins=linspace(xmin, xmax, 150), color='C0', label='disorder')
    axs[1].hist(order[feature_label], bins=linspace(xmin, xmax, 150), color='C1', label='order')
    axs[1].set_xlabel(f'Rate ({feature_label})')
    for i in range(2):
        axs[i].set_ylabel('Number of regions')
        axs[i].legend()
    plt.savefig(f'out/rates/{feature_label}.png')
    for i in range(2):
        axs[i].set_yscale('log')
    plt.savefig(f'out/rates/{feature_label}_log.png')
    plt.close()

# 3.2 Plot rate PCAs
pca = PCA(n_components=10)
idx = rate.index.get_level_values('disorder').array.astype(bool)

transform = pca.fit_transform(rate.to_numpy())
plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, color='C0', alpha=0.1, edgecolors='none')
plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, color='C1', alpha=0.1, edgecolors='none')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('unnorm')
legend = plt.legend(markerscale=2)
for lh in legend.legendHandles:
    lh.set_alpha(1)
plt.savefig(f'out/rates/pca_rate_unnorm.png')
plt.close()

x = (rate-rate.mean())/rate.std()
transform = pca.fit_transform(x.to_numpy())
plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, color='C0', alpha=0.1, edgecolors='none')
plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, color='C1', alpha=0.1, edgecolors='none')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('z-score')
legend = plt.legend(markerscale=2)
for lh in legend.legendHandles:
    lh.set_alpha(1)
plt.savefig(f'out/rates/pca_rate_z-score.png')
plt.close()

x = (rate-rate.min())/(rate.max()-rate.min())
transform = pca.fit_transform(x.to_numpy())
plt.scatter(transform[idx, 0], transform[idx, 1], label='disorder', s=5, color='C0', alpha=0.1, edgecolors='none')
plt.scatter(transform[~idx, 0], transform[~idx, 1], label='order', s=5, color='C1', alpha=0.1, edgecolors='none')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('min-max')
legend = plt.legend(markerscale=2)
for lh in legend.legendHandles:
    lh.set_alpha(1)
plt.savefig(f'out/rates/pca_rate_min-max.png')
plt.close()

# 3.3 Plot regions with extreme rates
for feature_label in rate.columns:
    if not os.path.exists(f'out/rates/{feature_label}/'):
        os.makedirs(f'out/rates/{feature_label}/')

    regions = set()
    ranked = df2[feature_label].sort_values(ascending=False).iteritems()
    while len(regions) < 20:
        (OGid, start, stop, _), _ = next(ranked)
        if (OGid, start, stop) in regions:
            continue
        regions.add((OGid, start, stop))

        # Load MSA, filter seqs, and re-order
        msa1 = load_msa(f'../insertion_trim/out/{OGid}.mfa')
        msa1 = {re.search(spid_regex, header).group(1): seq for header, seq in msa1}

        spids = region2spids[(OGid, start, stop)]
        tree = tree_template.shear(spids)
        order = {tip.name: i for i, tip in enumerate(tree.tips())}
        msa2 = [msa1[spid].upper()[start:stop] for spid in sorted(spids, key=lambda x: order[x])]
        fig = plot_msa(msa2, figsize=(8, 6), x_start=start)
        plt.savefig(f'out/rates/{feature_label}/{len(regions)-1}_{OGid}-{start}-{stop}.png', bbox_inches='tight')
        plt.close()

"""
DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/regions_30.tsv
../get_contrasts/get_contrasts.py
    ../get_contrasts/out/contrasts.tsv
"""