"""Plot qq-plots and histograms of the mixture models."""

import json
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pandas as pd
import scipy.stats as stats
from math import ceil
from numpy import linspace
from pymix.mixture import MixtureModel


def make_plots(pics_lt, dirpath, filename):
    # Load model
    with open(dirpath + '/' + filename) as file:
        dist_names, params, params_fix, weights, name = json.load(file)
    dists = [dists_dict[dist_name] for dist_name in dist_names]
    mixmod = MixtureModel(dists, params, params_fix, weights)

    # Calculate quantiles
    _, feature = dirpath.split('/')
    q_obs = pics_lt.loc[pics_lt[feature] != 0, feature].sort_values()
    q_dist = [mixmod.ppf((i - 0.5) / len(q_obs)) for i in range(1, len(q_obs) + 1)]

    for cutoff in [5, 2, 1, 0]:
        if cutoff == 0:
            sliver = slice(None)  # Slices do not behave nicely with zeroes
        else:
            idx = ceil(cutoff * len(q_obs) / 100)  # Ceil ensures the cutoff fraction is removed
            sliver = slice(idx, -idx)

        # Plot quantiles
        fig, ax = plt.subplots()
        ax.scatter(q_dist[sliver], q_obs[sliver], s=10, edgecolors='none')
        ax.plot(ax.get_xlim(), ax.get_xlim(), color='k', linewidth=1)  # Plot y=x line
        ax.set_xlabel('Theoretical Quantiles')
        ax.set_ylabel('Observed Quantiles')
        ax.set_title(f'{feature}: {name} with {cutoff}% Trimmed Tails')
        plt.savefig(f'out/{feature}/qq_{cutoff}_{name}.png')
        plt.close()

        # Plot histograms
        fig, ax = plt.subplots()
        _, bins, _ = ax.hist(q_obs[sliver], bins=50, density=True, color='white', linewidth=1, edgecolor='black')
        x = linspace(min(bins), max(bins), 1000)
        for i in range(len(mixmod.dists)):
            ax.plot(x, mixmod.pdf_comp(x, comp=i), label=mixmod.dists[i].name)
        ax.plot(x, mixmod.pdf(x), label='total')
        ax.set_xlabel(feature)
        ax.set_ylabel('Density')
        ax.set_title(f'{feature}: {name} with {cutoff}% Trimmed Tails')
        ax.legend()
        plt.savefig(f'out/{feature}/hist_{cutoff}_{name}.png')
        plt.close()


# Input variables
path = '../pic_calc/pics.tsv'
num_processes = int(os.environ['SLURM_NTASKS'])
lt = 32
dists_dict = {'laplace': stats.laplace,
              'norm': stats.norm}

if __name__ == '__main__':
    # Read data and filter
    pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
    pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
                   (~pics.index.get_level_values('ordered').array.astype(bool))]
    dirfile = [(pics_lt, dirpath, filename) for dirpath, _, filenames in os.walk('out')
               for filename in filenames if filename.endswith('.json')]

    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(make_plots, dirfile)

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/pics.tsv
./mixture_pic_calc.py
    ./out/*/*.json
"""