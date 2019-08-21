"""Execute the dists.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/dists.py ../feature_calc/', shell=True)

"""
NOTES
The shapes of the distributions for each feature appear largely the same between the original and shuffled sequences.
    In many cases, the mode is the same value, but its frequency is decreased.

DEPENDENCIES
../feature_plot_scripts/dists.py
../feature_calc/feature_calc.py
    ../feature_calc/features_*.tsv
"""