"""Plot outputs of fitting GO term logit models."""

import os

import matplotlib.pyplot as plt
import pandas as pd

pdidx = pd.IndexSlice

models = pd.read_table('out/models.tsv').set_index(['GOid', 'label'])
disorder = models.loc[pdidx[:, 'disorder'], :]
order = models.loc[pdidx[:, 'order'], :]
combined = models.loc[pdidx[:, 'combined'], :]

if not os.path.exists('out/'):
    os.mkdir('out/')

for data, data_label, color in [(disorder, 'disorder', 'C0'), (order, 'order', 'C1'), (combined, 'combined', 'C2')]:
    for feature_label in models.columns:
        plt.hist(data[feature_label], bins=50, label=data_label, color=color)
        plt.xlabel(feature_label.capitalize())
        plt.ylabel('Number of models')
        plt.legend()
        plt.savefig(f'out/hist_modelnum-{feature_label}_{data_label}.png')
        plt.close()

    plt.scatter(data['sensitivity'], data['specificity'], label=data_label, facecolor=color, edgecolor='none', alpha=0.25, s=12)
    plt.xlabel('Sensitivity')
    plt.ylabel('Specificity')
    leg = plt.legend(markerscale=1.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/scatter_spec-sens_{data_label}.png')
    plt.close()

"""
DEPENDENCIES
./fit_logits.py
    ./out/models.tsv
"""