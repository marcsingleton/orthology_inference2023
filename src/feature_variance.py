import os
import matplotlib.pyplot as plt
import pandas as pd
import re
from matplotlib import cm
from sys import argv


def merge(s_in, thresh):
    s_out = {}
    for idx in s_in.index:
        if s_in[idx] < thresh:
            s_out['other'] = s_out.get('other', 0) + s_in[idx]
        else:
            s_out[idx] = s_in[idx]
    return pd.Series(s_out)


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
    i = path[j0 + 1:j1]

    subsets = (df.loc[:, :], df.loc[T_idx, :], df.loc[F_idx, :])
    idxs = ('all', T_idx, F_idx)
    names = ('All', T_name, F_name)
    for subset, idx, name in zip(subsets, idxs, names):
        # Make output directories for feature sets
        cur_dir = f'{name}/'
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)  # Recursive folder creation

        # Get fractional variances and consolidate
        var = subset.var()
        vars_frac = merge(var / sum(var), 0.01).sort_values(ascending=False)

        # Create color map
        cmap = cm.get_cmap('GnBu')
        colors = [cmap(1 - k / len(vars_frac)) for k, _ in enumerate(vars_frac)]

        # Plot as pie chart
        plt.figure()
        plt.subplots_adjust(right=0.7)
        plt.pie(vars_frac, labels=vars_frac.index, colors=colors, labeldistance=None, autopct='%1.1f%%', pctdistance=1.15)
        plt.legend(loc='center left', bbox_to_anchor=(1.025, 0.5))
        plt.title(f'Fraction of Overall Variance by Feature:\n{name} Subsequences')
        plt.savefig(cur_dir + f'var_{i}.png')

        # Save var and var_frac to tsv
        df = pd.DataFrame({'var': var, 'var_frac': var / sum(var)}).sort_values(by='var', ascending=False)
        df.to_csv(cur_dir + f'var_{i}.tsv', sep='\t')
