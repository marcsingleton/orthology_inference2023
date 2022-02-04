"""Execute the tsne_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/tsne_separate.py ../feature_calc_shuffle/out/ 1 conserved Conserved Diverged', shell=True)

"""
NOTES
ALL
    At low cutoffs, the projections are nearly identical for both the conserved and diverged classes.
    At the highest cutoffs, the projections are still largely the same, but the clusters in the shuffled projections are more discrete while maintaining the overall fan-like structure seen in the unshuffled sequences.
        This difference is pronounced for the diverged sequences and subtle for the conserved sequences.

    As noted previously, the kappa-omega clusters are significantly reduced in both sequence classes at higher cutoffs.
MINMAX
    Like in the norm feature set, the shuffled projections are more uniform than are the unshuffled projections at high length cutoffs.
        This effect is more pronounced in the norm feature set as a whole and more pronounced in the diverged sequences within the minmax feature set.
    The shuffled projections also divide the central cluster into two distinct lobes.
NORM
    At high length cutoffs, the shuffled projections are nearly uniform whereas the unshuffled projections have regions of high and low density in the central cluster.
        The unshuffled sequences do not form discrete clusters outside of the usual kappa-omega clusters, however.
ZNORM
    At the highest length cutoffs, the patterns are similar to those observed for the other feature sets with the shuffled projections having a more uniform distribution of points within the central cluster.
    At moderate length cutoffs, the differences are more subtle.
        Instead of forming vertices around the edges of the central cluster, the shuffled sequences have distinct periphery clusters.
        However, within these periphery clusters, the distribution is uniform in contrast to the highly concentrated vertices.

DEPENDENCIES
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/out/features_*.tsv
"""