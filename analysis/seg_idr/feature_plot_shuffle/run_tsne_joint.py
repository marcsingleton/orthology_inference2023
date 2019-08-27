"""Execute the tsne_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc_shuffle/ ord dis Ordered Disordered', shell=True)

"""
NOTES

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_?.tsv
"""