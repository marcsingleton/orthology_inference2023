"""Execute the feature_variance.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/feature_variance.py ../sample_feats/ 3 ordered Ordered Disordered', shell=True)

"""
DEPENDENCIES
../../../src/feature_variance.py
../sample_feats/run_sample_feats.py
    ../feature_calc/features_*.tsv
"""