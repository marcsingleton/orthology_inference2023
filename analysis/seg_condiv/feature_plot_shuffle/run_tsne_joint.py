"""Execute the tsne_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc_shuffle/ con div Conserved Diverged', shell=True)

"""
NOTES
ALL
    No clear differences between shuffled and observed sequences for the all feature set.
        The shuffled sequences have a similar "constellation" of clusters observed in the unshuffled sequences.
MINMAX
    The patterns are similar to those in the norm feature set; however, the differences between the unshuffled and shuffled sequences appear somewhat more subtle.
NORM
    At low length cutoffs, the unshuffled and shuffled sequences are largely indistinguishable.
    At high length cutoffs, the conserved and diverged sequences more clearly segregate in a gradient in the central cluster.
        The central cluster is also more uniform and less punctate.
ZNORM
    Projections comparisons are the same as the minmax feature set relative to the norm feature set.

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_?.tsv
"""