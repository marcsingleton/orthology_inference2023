"""Execute the pca_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../feature_calc/ ord dis Ordered Disordered', shell=True)

"""
NOTES

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../feature_calc/feature_calc.py
    ../feature_calc/features_*.tsv
"""