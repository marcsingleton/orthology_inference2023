"""Fit mixture models to rate distributions."""

import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import json
import scipy.stats as stats
from pymix.estimate import cfes
from pymix.mixture import MixtureModel
from random import random


def get_rand_params(rand_maxes):
    rand_params = [{key: val * random() for key, val in rand_max.items()} for rand_max in rand_maxes]
    return rand_params


def fit_model(rates, feature, model):
    # Filter data and unpack variables
    data = rates.loc[rates[feature] != 0, feature].sort_values()
    name, dists = model

    # Make output directories for feature models
    cur_dir = f'out/{feature}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Get maxes for generation of random initials
    rand_maxes = []
    for dist in dists:
        cfe = cfes[dist.name]  # Get closed-form estimator
        n = max(1, int(len(data) * random()))
        sample_params = pd.DataFrame([cfe(np.random.choice(data, n)) for _ in range(500)])
        rand_max = (sample_params.mean() + num_std * sample_params.std()).to_dict()
        rand_maxes.append(rand_max)

    results = []
    excepts = []  # For numerical exceptions (nan, inf)
    while len(results) < num_init:
        mixmod = MixtureModel(dists, name=name, params=get_rand_params(rand_maxes))
        try:
            n, ll = mixmod.fit(data, max_iter=1000)
        except RuntimeError:  # Catch RuntimeErrors from failure to converge
            pass

        # Store output based on ll value
        if np.isnan(ll) or np.isinf(ll):
            excepts.append((n, ll, mixmod))
        else:
            results.append((n, ll, mixmod))

    # Store model with largest log-likelihood
    _, _, mixmod_max = max(results, key=lambda x: x[1])
    model = ([dist.name for dist in mixmod_max.dists],
             mixmod_max.params,
             mixmod_max.params_fix,
             mixmod_max.weights,
             name)
    with open(cur_dir + f'model_{name}.json', 'w') as file:
        json.dump(model, file)

    # Store EM metadata
    results.extend(excepts)
    ns, lls, cons = zip(*[(n, ll, mixmod.converged) for n, ll, mixmod in results])
    df = pd.DataFrame.from_dict({'n_iter': ns, 'll': lls, 'converged': cons})
    df.to_csv(cur_dir + f'meta_{name}.tsv', sep='\t', na_rep='nan')


# Input variables
path = '../mixture_filter/out/rates_filter.tsv'
lt = 32
num_processes = int(os.environ['SLURM_NTASKS'])
num_init = 10  # Number of initializations for each model
num_std = 10  # Number of standard deviations above mean for max of the random initials
models = [('lognorm3', [stats.lognorm, stats.lognorm, stats.lognorm]),
          ('lognorm2', [stats.lognorm, stats.lognorm]),
          ('glognorm2', [stats.gamma, stats.lognorm, stats.lognorm]),
          ('elognorm2', [stats.expon, stats.lognorm, stats.lognorm]),
          ('plognorm2', [stats.pareto, stats.lognorm, stats.lognorm]),
          ('llognorm2', [stats.levy, stats.lognorm, stats.lognorm])]

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Read data
    rates = pd.read_csv(path, sep='\t', index_col=0).dropna()
    featmods = [(rates, feature, model) for feature in rates for model in models]

    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(fit_model, featmods)

"""
DEPENDENCIES
../mixture_filter/mixture_filter.py
    ../mixture_filter/out/rates_filter.tsv
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""