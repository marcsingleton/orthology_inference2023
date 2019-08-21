"""Plot distributions of features for conserved and diverged subsequences jointly."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from shared import fsets
from sys import argv

feature_dir = argv[1]  # Feature directory must end in /
paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Split into separate conserved and divergent dataframes
    df = pd.read_csv(feature_dir + path, sep='\t', keep_default_na=False, index_col=[0, 1])

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = int(path[j0+1:j1])

    for key in ['all', 'minmax', 'znorm']:
        # Calculate features
        feat_all = fsets[key](df)

        # Make output directories for cutoffs and feature sets
        cur_dir = f'dists/{i}/{key}/'
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)  # Recursive folder creation

        # Plot distributions
        for feature in feat_all.columns:
            plt.figure()
            _, bins, _ = plt.hist(feat_all[feature], bins=25)
            plt.xlabel(feature)
            plt.ylabel('Count')
            plt.savefig(cur_dir + f'dist_{feature}_joint.png')
            plt.close()

            plt.figure()
            plt.hist(feat_all.loc['con'][feature], label='Conserved', bins=bins, alpha=0.5)
            plt.hist(feat_all.loc['div'][feature], label='Diverged', bins=bins, alpha=0.5)
            plt.xlabel(feature)
            plt.ylabel('Count')
            plt.legend()
            plt.savefig(cur_dir + f'dist_{feature}_separate.png')
            plt.close()
