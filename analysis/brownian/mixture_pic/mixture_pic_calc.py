"""Fit mixture models to PIC distributions with Laplace and Gaussian mixtures after removing zeroes."""

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


def fit_model(pics_lt, feature, model):
    # Filter data and unpack variables
    data = pics_lt.loc[pics_lt[feature] != 0, feature].sort_values()
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
path = '../pic_calc/pics.tsv'
num_processes = int(os.environ['SLURM_NTASKS'])
num_init = 10  # Number of initializations for each model
num_std = 20  # Number of standard deviations above mean for max of the random initials
models = [('norm4', [stats.norm, stats.norm, stats.norm, stats.norm]),
          ('norm3', [stats.norm, stats.norm, stats.norm]),
          ('norm2', [stats.norm, stats.norm]),
          ('laplace4', [stats.laplace, stats.laplace, stats.laplace, stats.laplace]),
          ('laplace3', [stats.laplace, stats.laplace, stats.laplace]),
          ('laplace2', [stats.laplace, stats.laplace])]
lt = 32

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Read data and filter
    pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
    pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
                   (~pics.index.get_level_values('ordered').array.astype(bool))]
    featmods = [(pics_lt, feature, model) for feature in pics_lt for model in models]

    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(fit_model, featmods)

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/pics.tsv
"""