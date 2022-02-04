"""Plot heatmaps of posterior probabilities for Laplace mixture models."""

import json
import matplotlib.pyplot as plt
import re
import os
import pandas as pd
import scipy.stats as stats
import matplotlib.cm as cm
from pymix.mixture import MixtureModel

# Input variables
path = '../pic_calc/out/pics.tsv'
lt = 32
dists_dict = {'laplace': stats.laplace}
height = 250

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
               (~pics.index.get_level_values('ordered').array.astype(bool))]

for feature in os.listdir('out'):
    model_paths = [x for x in os.listdir('out/' + feature) if re.match(r'model_laplace.\.json', x)]
    fig, axs = plt.subplots(len(model_paths), 1, figsize=(6, 3.5))
    fig.suptitle(f'{feature}:\nPosterior Probabilities of Mixture Model Components', y=0.95, size=10)
    fig.subplots_adjust(top=0.85)
    for i, model_path in enumerate(model_paths):
        # Load model, sorting components by scale and instantiating distributions from names
        with open('/'.join(['out', feature, model_path])) as file:
            model = json.load(file)
        model_params, name = model[:-1], model[-1]
        dist_names, params, params_fix, weights = [list(x) for x in zip(*sorted(zip(*model[:-1]), key=lambda y: y[1]['scale']))]
        dists = [dists_dict[dist_name] for dist_name in dist_names]
        mixmod = MixtureModel(dists, params, params_fix, weights, name)

        # Create heatmap
        data = pics_lt.loc[pics_lt[feature] != 0, feature].sort_values()
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

    plt.savefig('out/' + feature + '/prob_plot_laplace.png')
    plt.close()

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
./mixture_pic_calc.py
    ./out/*/*.json
"""