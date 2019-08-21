"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_joint.py ../feature_calc/', shell=True)

"""
NOTES
The shape of the projections appear less dependent on the length cutoff for this set of features than for the amino acid compositions.
    This is likely b/c these features are less dependent on single amino acids and there are more features, which "dilute" the effect of features that do dependent on single amino acids
As with the amino acid composition feature set, the effect of increasing the cutoff is largely diminished after a length of 3
    Although the shape is largely unchanged with increasing cutoffs, the spread appears less "sparse" as there are fewer overlapping points
ALL
    The projections for the fast_len set are largely uninformative since it is dominated by the net charge distribution of the conserved sequences
        The differences net charge distribution is in turn a result of differences in the length distribution; the net charge distribution, however, is not log transformed, resulting in its outsized influence on the PCA
    There are some clusters in the diverged subsequences at higher cutoffs
        Analyzing the data without the net charge feature would potentially reveal additional structures
LEN
    There are no substantial differences from the projections with the ALL feature set, which supports the previous conclusion that the net charge was the dominant component of the previous projection
        If the net charge distribution were log transformed, it is likely that the length would have a greater impact on the PCA
NET_CHARGE
    Beautifully pulls apart the conserved and diverged subsequences; examination of the PCs shows this is largely an effect of the length distribution
    The clusters noted in the diverged subsequences are now apparent along the x-axis; PC1 is largely pI, so these clusters are largely a reflection of the limited number of pIs for short sequences 
LEN+NET_CHARGE
    The subsequence classes still largely overlap
    PC1 is largely pI and PC2 is largely the ED ratio
        Clearly unnormalized variables can have outsized effects on the PCA
NORM
    3 large clusters are apparent; however, these are likely artifactual given kappa and omega return -1 for sequences where they are undefined (i.e. sequences w/o charged residues or prolines)
        See pca_ko.py for further discussion
ZNORM
    Like in the fast_norm feature set, the projections mostly overlap, except here the conserved region is much more pronounced
    The appropriateness of a z transform is dubious, as most of the feature distributions are highly asymmetric
Correctly normalizing the features, particularly if the distribution is not symmetric, will likely be an ongoing challenge
    Investigating clustering algorithms which do not employ distance-based metrics will ameliorate this
It is also worthwhile investigating other clustering algorithms even if they are distance-based (t-SNE, hierarchical) and other normalization strategies (min-max, mean, unit, correlation)
    A relevant alternative hypothesis is there are not distinct classes of sequences, and the subsequences exist as a continuum in feature space
        This appears to be the case for the conserved and diverged sequences; it is unclear if it applies to IDRs
Regardless future analyses should focus on IDRs; there is no hypothesis that there are classes of conserved or diverged subsequences

DEPENDENCIES
../../../src/feature_plot_scripts/pca_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""