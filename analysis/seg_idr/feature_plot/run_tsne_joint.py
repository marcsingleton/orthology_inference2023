"""Execute the tsne_joint.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc/ ord dis Ordered Disordered', shell=True)

"""
NOTES

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../feature_calc/feature_calc.py
    ../feature_calc/features_*.tsv
"""