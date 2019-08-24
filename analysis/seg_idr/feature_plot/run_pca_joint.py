"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_joint.py ../feature_calc/ ord dis Ordered Disordered', shell=True)

"""
NOTES

DEPENDENCIES
../../../src/feature_plot_scripts/pca_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""