"""Execute the tsne_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/tsne_separate.py ../feature_calc_shuffle/out/ 1 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    The projections are generally similar. However at high length cutoffs, the fan-like structure is less uniform and more clustered in the shuffled sequences.
        The subclusters are still arranged in a fan-like structure however.
        The subclusters are stronger in the disordered subsequences than in the order subsequences.
MINMAX
    The projections are largely identical at all length cutoffs.
        There is generally slightly more separation between the clusters in the shuffled sequences.
        Generally, it is difficult to determine if the shuffled projections are more uniform.
NORM
    The shuffled projections are more uniform and have fewer nodes than do the unshuffled projections.
        This effect is strongest at high length cutoffs but still noticeable at short length cutoffs.
ZNORM
    The shuffled projections are largely identical for both ordered and disordered subsequences.
        At high cutoffs and high perplexities, a band is visible in the shuffled projections that is absent from the observed projections.

DEPENDENCIES
../../../src/seg_scripts/plot/tsne_separate.py
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/out/features_*.tsv
"""