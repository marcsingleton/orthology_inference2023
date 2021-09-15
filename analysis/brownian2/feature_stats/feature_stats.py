"""Plot statistics associated with features."""

import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from numpy import linspace
from sklearn.decomposition import PCA


# Load features
features = pd.read_table('../get_features/out/features.tsv')
features.loc[features['kappa'] == -1, 'kappa'] = 1
features.loc[features['omega'] == -1, 'omega'] = 1

# Parse segments
rows = []
with open('../aucpred_filter/out/segments_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        for ppid in ppids.split(','):
            rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'disorder': disorder == 'True', 'ppid': ppid})
segments = pd.DataFrame(rows).merge(features, how='left', on=['OGid', 'start', 'stop', 'ppid'])
OGs = segments.groupby(['OGid', 'start', 'stop', 'disorder'])
mean = OGs.mean()

pdidx = pd.IndexSlice
disorder = mean.loc[pdidx[:, :, :, True], :]
order = mean.loc[pdidx[:, :, :, False], :]

if not os.path.exists('out/'):
    os.mkdir('out/')

# Feature histograms
for feature_label in mean.columns:
    fig, axs = plt.subplots(2, 1, sharex=True)
    xmin, xmax = mean[feature_label].min(), mean[feature_label].max()
    axs[0].hist(disorder[feature_label], bins=linspace(xmin, xmax, 75), color='C0', label='disorder')
    axs[1].hist(order[feature_label], bins=linspace(xmin, xmax, 75), color='C1', label='order')
    axs[1].set_xlabel(f'Mean {feature_label}')
    for i in range(2):
        axs[i].set_ylabel('Number of OGs')
        axs[i].legend()
    plt.savefig(f'out/hist_numOGs-{feature_label}.png')
    plt.close()


# Individual PCAs
pca = PCA(n_components=10)
plots = [(disorder, 'disorder', 'unnorm'),
         (order, 'order', 'unnorm'),
         ((disorder - disorder.mean()) / disorder.std(), 'disorder', 'z-score'),
         ((order - order.mean()) / order.std(), 'order', 'z-score'),
         ((disorder - disorder.min()) / (disorder.max() - disorder.min()), 'disorder', 'min-max'),
         ((order - order.min()) / (order.max() - order.min()), 'order', 'min-max')]
colors = ['#e15759', '#499894', '#59a14f', '#f1ce63', '#b07aa1', '#d37295', '#9d7660', '#bab0ac',
          '#ff9d9a', '#86bcb6', '#8cd17d', '#b6992d', '#d4a6c8', '#fabfd2', '#d7b5a6', '#79706e']

for data, data_label, norm_label in plots:
    color = 'C0' if data_label == 'disorder' else 'C1'

    # PCA without arrows
    transform = pca.fit_transform(data.to_numpy())
    plt.scatter(transform[:, 0], transform[:, 1], label=data_label, color=color, s=5, alpha=0.1, edgecolors='none')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(norm_label)
    legend = plt.legend(markerscale=2)
    for lh in legend.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/scatter_pca_{data_label}_{norm_label}.png')
    plt.close()

    # PCA with arrows
    plt.scatter(transform[:, 0], transform[:, 1], label=data_label, color=color, s=5, alpha=0.1, edgecolors='none')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(norm_label)

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    scale = (xmax + ymax - xmin - ymin) / 4
    projections = sorted(zip(data.columns, pca.components_[:2].transpose()), key=lambda x: x[1][0]**2 + x[1][1]**2, reverse=True)

    handles = []
    for i in range(len(colors)):
        feature_label, (x, y) = projections[i]
        handles.append(Line2D([], [], color=colors[i%len(colors)], label=feature_label))
        plt.annotate('', xy=(scale*x, scale*y), xytext=(0, 0), arrowprops={'headwidth': 5, 'headlength': 5, 'width': 0.5, 'color': colors[i%len(colors)]})
    plt.legend(handles=handles, fontsize=6, loc='center', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'out/scatter_pca-arrow_{data_label}_{norm_label}.png')
    plt.close()

    # Scree plot
    plt.bar(range(1, len(pca.explained_variance_ratio_)+1), pca.explained_variance_ratio_, label=data_label, color=color)
    plt.xlabel('Principal component')
    plt.ylabel('Explained variance ratio')
    plt.title(norm_label)
    plt.legend()
    plt.savefig(f'out/bar_{data_label}_{norm_label}.png')
    plt.close()

"""
DEPENDENCIES
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/segments_30.tsv
../get_features/get_features.py
    ../get_features/out/features.tsv
"""