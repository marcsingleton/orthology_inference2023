"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/pca_joint.py ../sample_feats/out/ 3 conserved Conserved Diverged', shell=True)

"""
NOTES
The projections are largely not characterized by the spokes found in the plots for the amino acid PCAs
    This is likely b/c these features are less dependent on single amino acids and there are more features, which "dilute" the effect of features that do dependent on single amino acids
ALL
    The projections for the fast_len set are largely uninformative since it is dominated by the net_charge and net_charge_P distributions of the conserved sequences
        The large variance of the charge distributions is a result of the differences in length between sequences
        The charge distributions correlate with length, but they are not log transformed, resulting in their outsize influence on the PCA
    As the length cutoff increases, SCD increases as a fraction of overall variance (as described in the run_feature_variance.py notes)
        This change in the relative contributions of each feature is likely why the shape of the projections is more parabolic at higher cutoffs
DROP_X
    Because features with an intrinsically larger scale can dominate the projection without necessarily revealing relevant substructure, dropping the most variable features was attempt to compensate for their contribution without imposing arbitrary normalizations
        This is likely not the most effective approach since we do not want to completely discard this information even if its contribution to the substructure is not commensurate with its fraction of the overall variance
    These plots are somewhat difficult to interpret because the removed feature is not apparent without referring to the model file
        An additional nuance is because the X most variable features were dropped, not every plot is necessarily removing the same features
            As noted previously, the relative contribution of SCD increases with the length cutoff
        This likely causes the appearance of the projections to dramatically vary between cutoffs
            Removing the most variable components does not necessarily decompress the projections if the most remaining most variable component has a small number of outliers that inflate its variance
                This is obvious in the complete feature set; a small number of highly charged sequences inflate the variance and the bounding box of the plot, thereby compressing the majority of the sequences into a small area
            The subtle interplay between the relative contributions to the variance of each feature as a function of the length cutoff combined with the stochasticity of the sampling process is likely the driving factor for the large variation in the shapes of the projections
    At higher cutoffs, the projections have a lattice structure, likely a result of the discrete nature of the RK and ED ratios
MINMAX
    This normalization almost exactly reproduces the projections in the norm feature set
        Because kappa and omega were excluded from normalization, the sequences with a -1 value for kappa or omega likely still dominate the variance, reproducing the same structure seen in the norm feature set
    Interestingly, at the highest cutoff, the three kappa-omega clusters are not visible in the minmax projection
        Likely at higher cutoffs, another feature with a large variance not present in the norm feature set exceeds the contributions from the relatively few sequences in the kappa-omega clusters
NORM
    3 large clusters are apparent; however, these are artifactual given kappa and omega return -1 for sequences where they are undefined (i.e. sequences w/o charged residues or prolines)
        The projections where points are colored by their kappa and omega values clearly demonstrate a direct correspondence between these classes and the clusters
ZNORM
    The projections mostly overlap, though some separation of the conserved and diverged sequences is visible
    The appropriateness of a z transform is dubious, however, as most of the feature distributions are highly asymmetric
Correctly normalizing the features, particularly if the distribution is not symmetric, will likely be an ongoing challenge
    Investigating clustering algorithms which do not employ distance-based metrics will ameliorate this
It is also worthwhile investigating other clustering algorithms even if they are distance-based (t-SNE, hierarchical) and correlation-based normalization strategies
    A relevant alternative hypothesis is there are not distinct classes of sequences, and the subsequences exist as a continuum in feature space
        This appears to be the case for the conserved and diverged sequences; it is unclear if it applies to IDRs
Regardless future analyses should focus on IDRs; there is no hypothesis that there are classes of conserved or diverged subsequences

DEPENDENCIES
../../../src/seg_scripts/plot/pca_joint.py
../sample_feats/run_sample_feats.py
    ../sample_feats/out/features_*.tsv
"""