"""Plot principal components of feature sets and distributions of those features for conserved and diverged subsequences jointly."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import re
from shared import n_components, fsets, labels
from sklearn.decomposition import PCA
from sys import argv

feature_dir = argv[1]  # Feature directory must end in /
paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Load data
    df = pd.read_csv(feature_dir + path, sep='\t', index_col=[0, 1])

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Get indices for plotting
    condiv_idx = df.index.get_level_values(0)
    con = 'con' == condiv_idx
    div = 'div' == condiv_idx
    k = df['kappa'] == -1
    o = df['omega'] == -1

    for key, fset in fsets.items():
        # Calculate features
        feat = fset(df)

        # Calculate PCA and transform data
        pca = PCA(n_components=n_components)
        trans = pca.fit_transform(feat.to_numpy())

        # Make output directories for feature sets
        cur_dir = f'pca_joint/{key}/'
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)  # Recursive folder creation

        # Plot PCAs
        # One panel
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.225)
        ax.scatter(trans[con, 0], trans[con, 1], s=2, alpha=0.1, label='Conserved')
        ax.scatter(trans[div, 0], trans[div, 1], s=2, alpha=0.1, label='Diverged')
        ax.set_title('PCA of Features\nin Conserved and Diverged Subsequences')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        leg = fig.legend(bbox_to_anchor=(0.5, 0.05), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'pca{i}_combined.png')
        plt.close()

        # Two panels
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
        ax1.scatter(trans[con, 0], trans[con, 1], s=2, alpha=0.1, label='Conserved', color='C0')
        ax2.scatter(trans[div, 0], trans[div, 1], s=2, alpha=0.1, label='Diverged', color='C1')
        ax1.set_title('PCA of Features\nin Conserved and Diverged Subsequences')
        ax2.set_xlabel('PC1')
        ax1.set_ylabel('PC2')
        ax2.set_ylabel('PC2')
        leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'pca{i}_separate.png')
        plt.close()

        # Color code by kappa and omega
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.225)
        for cond, label in labels.items():
            ko_idx = list(map(lambda x: x == cond, zip(k, o)))
            ax.scatter(trans[ko_idx, 0], trans[ko_idx, 1], s=2, alpha=0.1, label=label)
        ax.set_title('PCA of Features\nGrouped by Kappa and Omega Values')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'pca{i}_ko.png')
        plt.close()

        # Write model and model summary
        with open(cur_dir + f'model{i}_summary.txt', 'w') as file:
            file.write(f'Explained variance ratio of first {n_components} components\n')
            file.write(str(pca.explained_variance_ratio_) + '\n\n')
            file.write('Principle components in feature space (truncated to the largest 5)\n')
            for component in pca.components_:
                component = map(lambda x: round(x, 4), component)
                pairs = sorted(list(zip(feat.columns, component)), key=lambda x: abs(x[1]), reverse=True)[:5]
                file.write(str(pairs) + '\n')
        with open(cur_dir + f'model{i}.pickle', 'wb') as file:
            pickle.dump(pca, file)
