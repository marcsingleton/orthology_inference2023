"""Plot heatmaps of posterior probabilities."""

import json
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.stats as stats
import matplotlib.cm as cm
from pymix.mixture import MixtureModel

# Input variables
path = '../mixture_filter/out/rates_filter.tsv'
lt = 32
dists_dict = {'gamma': stats.gamma,
              'expon': stats.expon,
              'levy': stats.levy,
              'lognorm': stats.lognorm,
              'pareto': stats.pareto}
height = 250

# Read data
rates = pd.read_csv(path, sep='\t', index_col=0).dropna()

for feature in os.listdir('out'):
    model_paths = [x for x in os.listdir('out/' + feature) if x.endswith('.json')]
    fig, axs = plt.subplots(len(model_paths), 1, figsize=(6, 6))
    fig.subplots_adjust(top=0.875, bottom=0.075)
    fig.suptitle(f'{feature}:\nPosterior Probabilities of Mixture Model Components', y=0.95, size=10)
    for i, model_path in enumerate(model_paths):
        # Load model
        with open('/'.join(['out', feature, model_path])) as file:
            dist_names, params, params_fix, weights, name = json.load(file)
        dists = [dists_dict[dist_name] for dist_name in dist_names]
        mixmod = MixtureModel(dists, params, params_fix, weights, name)

        # Create heatmap
        data = rates.loc[rates[feature] != 0, feature].sort_values()
        expts = mixmod.posterior(data)

        # Plot heatmap
        ax = axs[i]
        ax.imshow(expts, vmin=0, vmax=1, aspect='auto', extent=[0, len(data), 0, len(mixmod.dists)])
        ax.set_aspect(0.035 * len(data))
        ax.set_title(mixmod.name, size=8)
        ax.tick_params(labelsize=7.5)
        ax.set_yticks([i + 0.5 for i in range(len(mixmod.dists))])
        ax.set_yticklabels([dist.name for dist in dists])
        if i != len(model_paths) - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel('Rank', size=8)

    # Plot colorbar
    cbar = fig.colorbar(cm.ScalarMappable(), ax=axs, fraction=0.025)
    cbar.set_label('Probability', size=8)
    cbar.ax.tick_params(labelsize=7.5)

    plt.savefig('out/' + feature + '/prob_plot.png')
    plt.close()

"""
DEPENDENCIES
../mixture_filter/mixture_filter.py
    ../mixture_filter/out/rates_filter.tsv
./mixture_rate_calc.py
    ./out/*/*.json
"""