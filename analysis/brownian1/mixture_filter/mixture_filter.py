"""Remove extremes from PIC distributions and plot removal statistics."""

import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
from pymix.mixture import MixtureModel

# Input variables
path = '../pic_calc/out/pics.tsv'
lt = 32
thresh = 0.95
dists_dict = {'laplace': stats.laplace}

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
               (~pics.index.get_level_values('ordered').array.astype(bool))]

# Read model selection into dictionary
model_paths = {}
with open('../mixture_pic/out/models.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.rstrip('\n').split('\t')
        model_paths[fields[0]] = fields[1]

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

rates = {}
counts = {}
fracs = {}
for feature in pics_lt:
    # Load model, sorting components by scale and instantiating distributions from names
    with open(f'../mixture_pic/out/{feature}/model_{model_paths[feature]}.json') as file:
            model = json.load(file)
    model_params, name = model[:-1], model[-1]
    dist_names, params, params_fix, weights = [list(x) for x in zip(*sorted(zip(*model[:-1]), key=lambda y: y[1]['scale']))]
    dists = [dists_dict[dist_name] for dist_name in dist_names]
    mixmod = MixtureModel(dists, params, params_fix, weights, name)

    # Remove extremes
    raw = pics_lt[feature]
    idx = np.flatnonzero(mixmod.posterior(raw)[-1] < thresh)
    clean = raw[idx]
    rates[feature] = (clean ** 2).groupby('block_id').mean()

    # Count extremes
    counts[feature] = len(raw) - len(clean)
    fracs[feature] = (len(raw) - len(clean)) / (len(raw) - (raw == 0).sum())  # Fraction of non-zero contrasts

    # Plot histograms of contrast counts in each block
    y = clean.groupby('block_id').count().value_counts()
    plt.bar(y.index, y)
    plt.title(f'{feature}: Contrast Counts in Blocks')
    plt.ylabel('Count')
    plt.xlabel('Number of Contrasts in Block')
    plt.xlim(0, 10)
    plt.savefig(f'out/bar_{feature}.png')
    plt.close()

# Plot counts and fractions of contrasts removed
plt.figure(figsize=(12, 4))
plt.bar(list(counts.keys()), list(counts.values()), width=0.5)
plt.title('Counts of Contrasts Removed')
plt.ylabel('Count')
plt.xlim(-1, 45)
plt.xticks(rotation=45, rotation_mode='anchor', horizontalalignment='right', fontsize=8)
plt.tick_params(axis='x', which='major', pad=1)
plt.tight_layout()
plt.savefig('out/counts.png')
plt.close()

plt.figure(figsize=(12, 4))
plt.bar(list(fracs.keys()), list(fracs.values()), width=0.5)
plt.title('Fraction of Non-zero Contrasts Removed')
plt.ylabel('Fraction')
plt.xlim(-1, 45)
plt.xticks(rotation=45, rotation_mode='anchor', horizontalalignment='right', fontsize=8)
plt.tick_params(axis='x', which='major', pad=1)
plt.tight_layout()
plt.savefig('out/fracs.png')
plt.close()

pd.DataFrame(rates).to_csv('out/rates_filter.tsv', sep='\t')

"""
DEPENDENCIES
../mixture_pic/mixture_pic_calc.py
    ../mixture_pic/out/models.tsv
    ../mixture_pic/out/*/model_*.json
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""