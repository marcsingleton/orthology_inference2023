"""Execute the feature_variance.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/feature_variance.py ../feature_calc/', shell=True)

"""
NOTES
As the cutoff increases, net_charge, net_charge_P, and SCD increase as a proportion of overall variance
    Unlike other features which are intrinsically length normalized or restricted to a certain range, these features can scale arbitrarily with length.

DEPENDENCIES
../../../src/feature_variance.py
../feature_calc/run_feature_calc.py
    ../feature_calc/features_*.tsv
"""