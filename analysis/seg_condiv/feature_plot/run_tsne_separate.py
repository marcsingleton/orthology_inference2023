"""Execute the tsne_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/tsne_separate.py ../sample_feats/out/ 3 conserved Conserved Diverged', shell=True)

"""
NOTES
Projecting the subsequences separately was successful in accentuating their substructures.
    This is consistent with the commentary on nested reductions in the PCA analysis; clusters were apparent in the combined analysis, and these clusters were a small fraction of the overall data.
ALL
    Charge dominates the structure of the conserved projection (as a result of their correlation with unscaled sequence length).
    At low cutoffs, the diverged subsequences form discrete clusters at all but the highest perplexities.
        It is unclear if these clusters are biologically relevant.
    At high cutoffs, these clusters merge into banded fan-like structures.
        They are similar, but slightly less structured, than the fans in the conserved subsequences.
MINMAX
    The patterns in structure by cutoff are largely the same as the norm feature set.
    The only significant difference is the central cluster in the conserved subsequences appears to have two lobes at all cutoffs (though it is more apparent at higher cutoffs).
        Two lobes are also visible in the central cluster of the diverged subsequence projections at the highest cutoffs.
NORM
    At low to medium cutoffs, there are three major clusters corresponding to values of kappa and omega like in the combined projections.
        In the conserved subsequences, these clusters are fairly uniform whereas in the diverged subsequences they display some substructure.
    At the highest cutoffs, undefined kappa and omega clusters are nearly absent in the conserved subsequences and dramatically reduced in size in the diverged subsequences.
        The central clusters are highly uniform, but both the conserved and diverged subsequences retain a small degree of structure.
ZNORM
    Its star-like topology is broadly similar to the combined projection.
        The conserved projections are more diffuse with little to no boundaries between the central and peripherary clusters; additionally, the points are more subtle.
        The diverged projections are more punctate along with the periphery and more clearly resolve the periphery clusters.
    As the length cutoffs increase, the periphery clusters decrease in size, disappearing at moderate values for the conserved subsequence projections.
        Additionally, the diverged projections are less punctate and more diffuse, similar to the conserved projections.            

DEPENDENCIES
../sample_feats/run_sample_feats.py
    ../sample_feats/out/features_*.tsv
"""