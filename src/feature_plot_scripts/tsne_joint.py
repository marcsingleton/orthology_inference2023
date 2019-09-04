"""Plot t-SNE of feature sets for two classes of subsequences jointly."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import re
from shared import fsets, labels
from sklearn.manifold import TSNE
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
        # Calculate features and initialize
        feat = fset(df)
        kl_divs = []

        # Make output directories for feature sets
        cur_dir = f'tsne_joint/{i}/{key}/'
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)  # Recursive folder creation

        for j in [5 * 2 ** x for x in range(8)]:
            # Calculate t-SNE and transform data
            tsne = TSNE(n_components=2, perplexity=j)
            trans = tsne.fit_transform(feat.to_numpy())
            kl_divs.append((j, tsne.kl_divergence_))

            # Plot t-SNEs
            # One panel
            fig, ax = plt.subplots()
            fig.subplots_adjust(bottom=0.225)
            ax.scatter(trans[T, 0], trans[T, 1], s=2, alpha=0.1, label=T_name)
            ax.scatter(trans[F, 0], trans[F, 1], s=2, alpha=0.1, label=F_name)
            ax.set_title(f't-SNE of Features\nin {T_name} and {F_name} Subsequences')
            leg = fig.legend(bbox_to_anchor=(0.5, 0.05), loc='lower center', ncol=2, markerscale=2.5)
            for lh in leg.legendHandles:
                lh.set_alpha(1)
            plt.savefig(cur_dir + f'tsne{j}_combined.png')
            plt.close()

            # Two panels
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
            ax1.scatter(trans[T, 0], trans[T, 1], s=2, alpha=0.1, label=T_name, color='C0')
            ax2.scatter(trans[F, 0], trans[F, 1], s=2, alpha=0.1, label=F_name, color='C1')
            ax1.set_title(f't-SNE of Features\nin {T_name} and {F_name} Subsequences')
            leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
            for lh in leg.legendHandles:
                lh.set_alpha(1)
            plt.savefig(cur_dir + f'tsne{j}_separate.png')
            plt.close()

            # Color code by kappa and omega
            fig, ax = plt.subplots()
            fig.subplots_adjust(bottom=0.225)
            for cond, label in labels.items():
                ko_idx = list(map(lambda x: x == cond, zip(k, o)))
                ax.scatter(trans[ko_idx, 0], trans[ko_idx, 1], s=2, alpha=0.1, label=label)
            ax.set_title('t-SNE of Features\nGrouped by Kappa and Omega Values')
            leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
            for lh in leg.legendHandles:
                lh.set_alpha(1)
            plt.savefig(cur_dir + f'tsne{j}_ko.png')
            plt.close()

            # Write model
            with open(cur_dir + f'model{j}.pickle', 'wb') as file:
                pickle.dump(tsne, file)

        # Write model summaries
        with open(cur_dir + 'model_summary.txt', 'w') as file:
            file.write('#perplexity\tKL_divergence\n')
            for perp, kl_div in kl_divs:
                file.write(f'{perp}\t{kl_div}\n')
