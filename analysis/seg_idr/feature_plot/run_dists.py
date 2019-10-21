"""Execute the dists.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/dists.py ../sample_feats/ 3 ordered Ordered Disordered', shell=True)

"""
DEPENDENCIES
../../../src/feature_plot_scripts/dists.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""