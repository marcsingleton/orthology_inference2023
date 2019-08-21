"""Execute the tsne_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_separate.py ../feature_calc_shuffle/', shell=True)

"""
NOTES

DEPENDENCIES
../feature_plot_scripts/tsne_separate.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_4.tsv
"""