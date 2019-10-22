"""Execute the pca_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../feature_calc_shuffle/ 1 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    Generally the ordered projections are highly similar and fan-shaped, though the shuffled subsequence projections are increasingly wider at larger length cutoffs.
        The shuffled projections are larger, due to the sequences having a smaller range in the shuffled projections.
            The increase in size displays the discrete nature of the net charge feature, creating "stripes" in the projection.
        The observed ordered subsequences do not have a parabolic shape except in one case, interestingly at the highest cutoff.
            It appears the ordered subsequences generally have less extreme feature values than the disordered subsequences.
            The one example of the observed ordered subsequence with a parabolic projection suggests that the projections are somewhat stochastic and change dramatically if a few outlier subsequences are sampled.
    The shuffled disordered subsequences are uniformly fan-shaped unlike the parabolic shapes of their observed counterparts.
        At higher length cutoffs, the distributions grow wider as observed with the ordered subsequences.
        Unlike the ordered subsequences, the projections of the disordered subsequences are more irregular at higher cutoffs.
MINMAX
    In both the observed and shuffled sequences, kappa-omega clusters are visible, but the shuffled subsequences generally have additional substructure such as partially resolved subclusters or bands.
        However, at some cutoffs, these are also visible in the observed sequences as well.
    At low cutoffs, these distinctions are largely subtle and inconsistent. However at high cutoffs, the shuffled subsequences clearly have a banded structure that is absent in the observed subsequences.
        Generally, both projections appear as "smeared" clusters, however.
NORM
    All projections are nearly identical except at the highest cutoffs where the shuffled projections are less angular and more uniformly Gaussian.
ZNORM
    At low cutoffs, the projections are largely identical.
    At high cutoffs, the shuffled projections are more uniform and Gaussian whereas the observed projections are angular and compressed.

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""