"""Execute the pca_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/pca_separate.py ../feature_calc_shuffle/out/ 1 conserved Conserved Diverged', shell=True)

"""
NOTES
ALL
    Differences between shuffled and observed sequences similar to those noted in run_pca_joint.py for the all feature set.
MINMAX
    The projections largely follow the same patterns observed for the norm feature set.
    The kappa-omega clusters tend to disappear at lower cutoffs than in the norm feature set projections.
        The sequences are still present of course, but their contribution to the variance is not sufficient to ensure their separation.
        This shows that although kappa-omega clusters are not visible in the projection, they are still potentially present in significant numbers.
            The best test for enrichment of kappa-omega clusters over expectation is by counting and performing a goodness-of-fit test.
NORM
    The separate projections for the conserved and diverged classes highlights in the differences with regards to the presence of the kappa and omega clusters.
        The kappa and omega clusters disappear by a length cutoff of 16 in the shuffled conserved sequences.
            The clusters are still present in the original sequences until a cutoff of 32.
        In contrast, the clusters do not disppear until a length cutoff of 32 in the shuffled diverged sequences.
        In both cases, the kappa and omega clusters are denser (or more enriched) in the observed sequences.
            It is unclear if the kappa and omega clusters are more enriched in diverged than conserved sequences relative to their shuffled counterparts.
                Because the shuffled diverged sequences cluster by kappa and omega values at a higher cutoff, the amino acid pool is likely more biased towards producing sequences with no prolines or charged residues.
                It is unclear if the increased number of sequences in kappa and omega clusters outstrips their inherent propensity for the kappa and omega clusters.
                Directly counting the sequences in all categories should more clearly distinguish these possibilities.
ZNORM
    Interestingly, even at low cutoffs, the conserved unshuffled sequences are largely uniform in their projections whereas the conserved shuffled sequences have bands of increased intensity.
        These largely disappear after the first few cutoffs.
    The shuffled diverged sequences are intensely banded compared with their unshuffled counterparts, which are characterized by spokes.
        Again these differences largely disppear at moderate cutoffs.
        However, the unshuffled sequences still have remnants of spokes whereas the shuffled sequences are nearly uniform.
        These observations are consistent with the conclusion that shuffled sequences are more "average" and observed sequences have a higher proportion of outliers.

DEPENDENCIES
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/out/features_*.tsv
"""