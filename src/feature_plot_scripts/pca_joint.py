"""Plot principal components of feature sets for two classes of subsequences jointly."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import re
from shared import n_components, fsets, labels
from sklearn.decomposition import PCA
from sys import argv

# Input variables
feature_dir = argv[1]  # Feature directory must end in /
T_idx = argv[2]  # Index of True class in sentence case
F_idx = argv[3]  # Index of False class in sentence case
T_name = argv[4]  # Name of True class in sentence case
F_name = argv[5]  # Name of False class in sentence case

paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Load data
    df = pd.read_csv(feature_dir + path, sep='\t', index_col=[0, 1])

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Get indices for plotting
    TF_idx = df.index.get_level_values(0)
    T = T_idx == TF_idx
    F = F_idx == TF_idx
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
        ax.scatter(trans[T, 0], trans[T, 1], s=2, alpha=0.1, label=T_name)
        ax.scatter(trans[F, 0], trans[F, 1], s=2, alpha=0.1, label=F_name)
        ax.set_title(f'PCA of Features\nin {T_name} and {F_name} Subsequences')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        leg = fig.legend(bbox_to_anchor=(0.5, 0.05), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'pca{i}_combined.png')
        plt.close()

        # Two panels
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
        ax1.scatter(trans[T, 0], trans[T, 1], s=2, alpha=0.1, label=T_name, color='C0')
        ax2.scatter(trans[F, 0], trans[F, 1], s=2, alpha=0.1, label=F_name, color='C1')
        ax1.set_title(f'PCA of Features\nin {T_name} and {F_name} Subsequences')
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
