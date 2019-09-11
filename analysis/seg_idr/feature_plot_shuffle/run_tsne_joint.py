"""Execute the tsne_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc_shuffle/ ord dis Ordered Disordered', shell=True)

"""
NOTES
ALL
    The projections are largely identical at all cutoffs.
        At the highest cutoffs, however, the shuffled sequences more strongly form "triangular" subclusters arranged in a fan-like pattern.
        In contrast, the lower cutoff projections form a more uniform spread within the fan-like pattern.
MAXMIN
    At low cutoffs, the projections are nearly identical.
    At high cutoffs (16 and 32), the projections have dramatic differences.
        In the observed sequences, the ordered and disordered sequences form a gradient and another gap begins to divide this gradient, forming 4 loosely connected subclusters.
            These 4 subclusters touch at their exteriors, creating a gap of no density in the center.
        In the shuffled sequences, this gradient and bisecting gap is still present. However, the 4 subclusters touch at their center, creating 4 "stumps" that project from a connected region in the center.
        The contrast is largest at a cutoff of 16. At a cutoff of 32, shuffled subclusters clearly separate.
            The observed subsequences, particularly the disordered subsequences, do not fully separate. However, the 4 cluster structure is still more obvious than in the 16 cutoff projections.
    What these subclusters are is unclear. They are not kappa-omega clusters.
        It is potentially iso_point, as the distributions for iso_point are bimodal for both ordered and disordered subsequences.
NORM
    At low cutoffs, the projections are identical.
    At high cutoffs, the shuffled subsequences are more uniform and display few to no variations in density. The overall structure from the observed subsequences is preserved, however.
ZNORM
    The projections are nearly indistinguishable at high and low cutoffs.
        At high cutoffs, the shuffled projections are potentially more uniform as a whole, but they also appear to contain a few more bands and nodes than the observed projections.

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_?.tsv
"""