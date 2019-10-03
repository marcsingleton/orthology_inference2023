"""Plot distributions of features for two classes of subsequences jointly and separately."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from sys import argv

# Input variables
feature_dir = argv[1]  # Feature directory must end in /
T_idx = argv[2]  # Index of True class
F_idx = argv[3]  # Index of False class
T_name = argv[4]  # Name of True class
F_name = argv[5]  # Name of False class

paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Load data
    df = pd.read_csv(feature_dir + path, sep='\t', index_col=[0, 1])

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = int(path[j0+1:j1])

    # Make output directories for cutoffs and feature sets
    cur_dir = f'dists/{i}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Plot distributions
    for feature in df.columns:
        # Plot classes as one series
        plt.figure()
        _, bins, _ = plt.hist(df[feature], bins=25)
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.savefig(cur_dir + f'{feature}_joint.png')
        plt.close()

        # Plot classes as two series
        plt.figure()
        plt.hist(df.loc[T_idx][feature], label=T_name, bins=bins, alpha=0.5)
        plt.hist(df.loc[F_idx][feature], label=F_name, bins=bins, alpha=0.5)
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_separate.png')
        plt.close()

        # Plot classes as separate graphs
        plt.figure()
        plt.hist(df.loc[T_idx][feature], label=T_name, bins=25, color='C0')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{T_idx}.png')
        plt.close()

        plt.figure()
        plt.hist(df.loc[F_idx][feature], label=F_name, bins=25, color='C1')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{F_idx}.png')
        plt.close()
