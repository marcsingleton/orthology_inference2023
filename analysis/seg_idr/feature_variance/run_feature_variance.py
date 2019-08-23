"""Execute the feature_variance.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/feature_variance.py ../feature_calc/', shell=True)

"""
DEPENDENCIES
../../../src/feature_variance.py
../feature_calc/run_feature_calc.py
    ../feature_calc/features_*.tsv
"""