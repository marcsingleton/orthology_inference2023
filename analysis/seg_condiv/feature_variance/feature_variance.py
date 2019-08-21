import os
import matplotlib.pyplot as plt
import pandas as pd
import re
from matplotlib import cm


def merge(s_in, thresh):
    s_out = {}
    for idx in s_in.index:
        if s_in[idx] < thresh:
            s_out['other'] = s_out.get('other', 0) + s_in[idx]
        else:
            s_out[idx] = s_in[idx]
    return pd.Series(s_out)


feature_dir = '../feature_calc/'
paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
for path in paths:
    # Load data
    df = pd.read_csv(feature_dir + path, sep='\t', index_col=[0, 1])

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Get fractional variances and consolidate
    var = df.var()
    vars_frac = merge(var / sum(var), 0.01).sort_values(ascending=False)

    # Create color map
    cmap = cm.get_cmap('GnBu')
    colors = [cmap(1 - i / len(vars_frac)) for i, _ in enumerate(vars_frac)]

    # Plot as pie chart
    plt.figure()
    plt.subplots_adjust(right=0.7)
    plt.pie(vars_frac, labels=vars_frac.index, colors=colors, labeldistance=None, autopct='%1.1f%%', pctdistance=1.15)
    plt.legend(loc='center left', bbox_to_anchor=(1.025, 0.5))
    plt.title('Fraction of Overall Variance by Feature')
    plt.savefig(f'var_{i}.png')
