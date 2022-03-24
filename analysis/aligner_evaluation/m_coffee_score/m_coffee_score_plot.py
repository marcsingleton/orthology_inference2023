"""Create plots from benchmark scores."""

import os

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA


def pca_plots(df, metric):
    S = df.drop(['clustalo2', 'clustalo3'], level=2)[metric].unstack(2)
    pca = PCA(n_components=5)
    pca.fit(S.to_numpy())

    # 1.1 PC1 and PC2
    plt.figure(dpi=400)
    for i, (ref, group) in enumerate(S.groupby('ref')):
        T = pca.transform(group)
        plt.scatter(T[:, 0], T[:, 1], label=ref, s=20, alpha=0.75, color=f'C{i}', edgecolors='none')
    plt.legend()
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(metric[:-1] + ' score')
    plt.savefig(f'out/scatter_pca12_{metric}.png')
    plt.close()

    # 1.2 PC2 and PC3
    plt.figure(dpi=400)
    for i, (ref, group) in enumerate(S.groupby('ref')):
        T = pca.transform(group)
        plt.scatter(T[:, 1], T[:, 2], label=ref, s=20, alpha=0.75, color=f'C{i}', edgecolors='none')
    for i, col in enumerate(S.columns):
        x = pca.components_[1, i] / 5
        y = pca.components_[2, i] / 5
        if x > 0:
            ha = 'left'
            x_offset = 0.006
        elif x == 0:
            ha = 'center'
            x_offset = 0
        else:
            ha = 'right'
            x_offset = -0.006
        plt.arrow(0, 0, x, y, width=0.0001, head_width=0.003)
        plt.text(x + x_offset, y, col, fontsize='small', ha=ha)
    plt.legend()
    plt.xlabel('PC2')
    plt.ylabel('PC3')
    plt.title(metric[:-1] + ' score')
    plt.savefig(f'out/scatter_pca23_{metric}.png')
    plt.close()


def scatter_matrix(df, metric):
    S = df.drop(['clustalo2', 'clustalo3'], level=2)[metric].unstack(2)
    fig, axs = plt.subplots(len(S.columns)-1, len(S.columns)-1, sharex=True, sharey=True)
    axs[0, 0].set_xlim(-0.1, 1.1)
    axs[0, 0].set_xticks([0, 1])
    axs[0, 0].set_ylim(-0.1, 1.1)
    axs[0, 0].set_yticks([0, 1])
    for i in range(len(S.columns)-1):
        for j in range(len(S.columns)-1):
            if j > i:
                axs[i, j].axis('off')
                continue
            if j == 0:
                axs[i, j].set_ylabel(S.columns[i+1], size='small')
            if i == len(S.columns)-2:
                axs[i, j].set_xlabel(S.columns[j], size='small')
            axs[i, j].scatter(S[S.columns[i+1]], S[S.columns[j]], s=5, alpha=0.5, edgecolors='none')
    fig.suptitle(metric[:-1] + ' score')
    fig.savefig(f'out/scatter_matrix_{metric}.png')
    plt.close()


df1 = pd.read_table('out/scores.tsv').set_index(['ref', 'file_id', 'aligner'])
df2 = pd.read_table('../bm_score/out/scores.tsv').set_index(['ref', 'file_id', 'aligner'])
df = pd.concat([df1, df2])
df['TC1'] = df['test_columns1'] / df['ref_columns1']
df['TC2'] = df['test_columns2'] / df['ref_columns2']
df['SP1'] = df['test_pairs1'] / df['ref_pairs1']
df['SP2'] = df['test_pairs2'] / df['ref_pairs2']

if not os.path.exists('out/'):
    os.mkdir('out/')

summary = df.groupby(['ref', 'aligner'])[['TC1', 'TC2', 'SP1', 'SP2']].agg(['mean', 'median', 'std', 'min', 'max'])
summary.to_csv('out/summary.tsv', sep='\t')

scatter_matrix(df, 'SP2')
scatter_matrix(df, 'TC2')
pca_plots(df, 'SP2')
pca_plots(df, 'TC2')

"""
DEPENDENCIES
../bm_score/bm_score_calc.py
    ../bm_score/out/scores.tsv
./m_coffee_score_calc.py
    ./out/scores.tsv
"""