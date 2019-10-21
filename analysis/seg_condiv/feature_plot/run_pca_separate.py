"""Execute the pca_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../sample_feats/ 3 conserved Conserved Diverged', shell=True)

"""
NOTES
The principal components for the conserved and diverged subsequences were analyzed separately to reveal any substructures in either that was obscured by analyzing them together
    The projections were largely the same
Nested dimensionality reductions are best for revealing substructure when classes are apparent in the primary analysis and the classes are sparse relative to the data set
    This data set fails both requirements; there are two equally-sized classes that largely overlap
    In contrast, many secondary t-SNE analyses in single cell data sets reveal sub cell types because:
        1. Clusters (cell types) were apparent in the primary reduction
        2. The number of points in each cluster is small relative to the entire data set
ALL
    Banding is particularly strong in the diverged sequences with smaller length cutoffs
    Since the first principal component is largely a combination of net_charge, net_charge_P, and iso_point, this is likely an effect of the discrete nature of charged residues
DROP_X
    No distinct clusters are visible; furthermore, the projections are difficult to interpret for many of the same reasons as discussed in run_pca_joint.py
MINMAX/NORM
    The separate projections are largely the same as the joint projections
    In both the minmax and norm feature sets, there are kappa-omega clusters in the diverged sequences but no clusters in the conserved sequences at the highest length cutoff
        Does this represent an enrichment of these sequences, or do the length and amino acid distributions of the diverged sequences have an intrinsic propensity for subsequences with undefined values for kappa and omega?
            Comparing these projections to those of the shuffled sequences will clarify
    In the minmax feature set only, it appears there are two distinct but partially overlapping populations of conserved sequences
ZNORM
    Spokes are visible in the diverged sequence projections at lower length cutoffs
    At larger length cutoffs, the spokes largely disappear
        This is a pattern common across all the feature sets and projections
            At lower length cutoffs, the projections (particularly those for the diverged sequences) are more asymmetrical or possess bands
                These structures are likely the result of large numbers of short and nearly identical sequences
            At higher length cutoffs, the projections are more symmetrical and diffuse
                As the number of possible sequences increase, the sequences are more uniformly distributed in sequence space

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""