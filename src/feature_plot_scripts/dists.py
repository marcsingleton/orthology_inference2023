"""Plot distributions of features for two types of subsequences jointly and separately."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from sys import argv

# Input variables
feature_dir = argv[1]  # Feature directory must end in /
num_idx = int(argv[2])  # Number of index columns
type_name = argv[3]  # Name of column denoting segment type
T_name = argv[4]  # Name of True type in sentence case
F_name = argv[5]  # Name of False type in sentence case

paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Read data
    df = pd.read_csv(feature_dir + path, sep='\t', index_col=list(range(num_idx)))

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = int(path[j0+1:j1])

    # Get indices for plotting
    T_idx = df.index.get_level_values(type_name).array.astype(bool)
    F_idx = ~T_idx

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
        plt.hist(df.loc[T_idx, feature], label=T_name, bins=bins, alpha=0.5)
        plt.hist(df.loc[F_idx, feature], label=F_name, bins=bins, alpha=0.5)
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_separate.png')
        plt.close()

        # Plot classes as separate graphs
        plt.figure()
        plt.hist(df.loc[T_idx, feature], label=T_name, bins=25, color='C0')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{T_name}.png')
        plt.close()

        plt.figure()
        plt.hist(df.loc[F_idx, feature], label=F_name, bins=25, color='C1')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{F_name}.png')
        plt.close()
