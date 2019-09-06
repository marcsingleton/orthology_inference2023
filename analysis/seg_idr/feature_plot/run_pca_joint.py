"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_joint.py ../feature_calc/ ord dis Ordered Disordered', shell=True)

"""
NOTES
ALL
    Across all length cutoffs, the projections have the parabolic shape observed for the longer cutoffs of the condiv subsequences.
MINMAX
    At cutoffs >= 16, the kappa-omega clusters are not observed.
    The subsequences do not separate into classes as observed for the norm feature set.
NORM
    At all but the high length cutoff, the projections have the characteristic kappa-omega clusters.
    At the highest length cutoff, the subsequences clearly separate into ordered and disordered classes.
        This is potentially an effect of the decreased number of kappa-omega sequences rather than a sudden increase in the differences between the ordered and disordered classes.
ZNORM
    The lower length cutoff projections have fewer spokes than their counterparts in the condiv subsequences.
        This is potentially a result of the average longer lengths of the subsequences in this set reducing the number of mono, di, and tri peptides.
    At length cutoffs >= 16, the two classes begin to clearly separate.

DEPENDENCIES
../../../src/feature_plot_scripts/pca_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""