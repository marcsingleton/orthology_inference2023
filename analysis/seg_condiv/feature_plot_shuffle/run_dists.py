"""Execute the dists.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/dists.py ../feature_calc_shuffle/ 1 conserved Conserved Diverged', shell=True)

"""
DEPENDENCIES
../../../src/feature_plot_scripts/dists.py
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""