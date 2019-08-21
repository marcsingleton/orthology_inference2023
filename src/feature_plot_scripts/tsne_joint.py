"""Plot t-SNE of feature sets for conserved and diverged subsequences jointly."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import re
from shared import fsets, labels
from sklearn.manifold import TSNE
from sys import argv

# Load data
feature_dir = argv[1]  # Feature directory must end in /
paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
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

            # Plot PCAs
            # One panel
            fig, ax = plt.subplots()
            fig.subplots_adjust(bottom=0.225)
            ax.scatter(trans[con, 0], trans[con, 1], s=2, alpha=0.1, label='Conserved')
            ax.scatter(trans[div, 0], trans[div, 1], s=2, alpha=0.1, label='Diverged')
            ax.set_title('t-SNE of Features\nin Conserved and Diverged Subsequences')
            leg = fig.legend(bbox_to_anchor=(0.5, 0.05), loc='lower center', ncol=2,markerscale=2.5)
            for lh in leg.legendHandles:
                lh.set_alpha(1)
            plt.savefig(cur_dir + f'tsne{j}_combined.png')
            plt.close()

            # Two panels
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
            ax1.scatter(trans[con, 0], trans[con, 1], s=2, alpha=0.1, label='Conserved', color='C0')
            ax2.scatter(trans[div, 0], trans[div, 1], s=2, alpha=0.1, label='Diverged', color='C1')
            ax1.set_title('t-SNE of Features\nin Conserved and Diverged Subsequences')
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
