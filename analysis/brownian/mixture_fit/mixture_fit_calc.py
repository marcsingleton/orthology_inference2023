"""Fit mixture models to rate distributions."""

import os
import pandas as pd
import json
import scipy.stats as stats
from mixture import MixtureModel

# Input variables
path = '../pic_calc/pics.tsv'
lt = 32
n_iter = 5
models = [('lognorm3', [stats.lognorm, stats.lognorm, stats.lognorm], [['s', 'scale'], ['s', 'scale'], ['s', 'scale']]),
          ('lognorm2', [stats.lognorm, stats.lognorm], [['s', 'scale'], ['s', 'scale']]),
          ('glognorm2', [stats.gamma, stats.lognorm, stats.lognorm], [['a', 'scale'], ['s', 'scale'], ['s', 'scale']]),
          ('elognorm2', [stats.expon, stats.lognorm, stats.lognorm], [['scale'], ['s', 'scale'], ['s', 'scale']]),
          ('plognorm2', [stats.pareto, stats.lognorm, stats.lognorm], [['b', 'scale'], ['s', 'scale'], ['s', 'scale']]),
          ('llognorm2', [stats.levy, stats.lognorm, stats.lognorm], [['scale'], ['s', 'scale'], ['s', 'scale']])]

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) & (~pics.index.get_level_values('ordered').array.astype(bool))]
rates = (pics_lt ** 2).groupby('block_id').mean()

for feature in rates:
    data = rates.loc[rates[feature] != 0, feature].sort_values()

    # Make output directories for feature models
    cur_dir = f'out/{feature}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Fit models
    for name, dists, params in models:
        print(feature, name, sep=': ')

        results = []
        for i in range(n_iter):
            mixmod = MixtureModel(dists, params)
            try:
                mixmod.fit(data)
                results.append(mixmod)
            except RuntimeError:  # Catch RuntimeErrors from failure to converge
                pass

        # Store model with largest log-likelihood
        max_ll = max(results, key=lambda x: x.ll)
        model = ([dist.name for dist in max_ll.dists], max_ll.params_fix, max_ll.params, max_ll.weights)
        with open(cur_dir + f'model_{name}.json', 'w') as file:
            json.dump(model, file)

        # Store EM metadata
        lls, ns, cons = zip(*[(mixmod.ll, mixmod.n_iter, mixmod.converged) for mixmod in results])
        df = pd.DataFrame.from_dict({'ll': lls, 'n_iter': ns, 'converged': cons})
        df.to_csv(cur_dir + f'meta_{name}.tsv', sep='\t')
