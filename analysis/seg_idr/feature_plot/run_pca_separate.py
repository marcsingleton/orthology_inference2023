"""Execute the pca_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../sample_feats/ 3 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    The ordered sequences generally have a fan shape while the disordered sequences have a parabolic shape.
        The parabolic shape is likely an artifact of the SCD feature, and largely driven by a few outliers in the disordered class which inflate the variance.
        Often projections excluding the most variable feature have a fan shape. In these cases we are likely viewing largely the same distribution along a different axis.
MINMAX
    In both classes, the kappa-omega clusters are visible at low to moderate cutoffs, but they are more "smeared" as compared to the norm feature set.
        They largely disappear at higher cutoffs.
NORM
    The kappa-omega clusters are visible in both classes. While smaller at larger cutoffs, they do not disappear completely.
ZNORM
    As noted in run_pca_joint.py, the lower length cutoffs have fewer spokes compared to the condiv projections.
    Both classes largely have no structure and the projections do not vary greatly with length cutoff.

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""