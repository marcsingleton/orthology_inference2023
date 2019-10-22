"""Execute the tsne_joint.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../sample_feats/ 3 conserved Conserved Diverged', shell=True)

"""
NOTES
The perplexity has a significant impact on the topology of the projection
    At the lowest perplexity, the projections are always "balls" of uniform density; though, conserved and diverged subsequences form distinct regions
    Significant structures do not tend to emerge until the perplexity is 20
    The projections at the highest perplexity are usually highly structured. None appear to have an entirely uniform "ball" topology, indicating the perplexity is not set too high
        At higher length cutoffs, the substructure is smoother and more subtle, likely as a result of the increased lengths allowing the sequences to cover more feature space
Generally substructures are apparent in the clusters, even if the clusters themselves are largely uniform
    However, it is not clear if the clusters and bands are an artifact of the t-SNE projection and feature set, or if they represent biologically distinct subsequences
ALL
    The conserved subsequences form banded fans which converge to a point; these are likely clustered by some net charge distribution
    The diverged subsequences form a series of clusters along the periphery
        Many are small, on the order of tens of points, and poorly resolved from adjacent clusters
        Some are larger, but these often have substructure as well suggesting even broadly similar sequences are composed of heterogeneous subpopulations
        Though the size and poor resolution of the clusters are poor, they are largely more defined than simulations with multidimensional Gaussians
    The conserved and diverged subsequences are not completely resolved. In particular, though they are largely composed of diverged subsequences, the periphery clusters also contain some conserved subsequences
    The periphery clusters are highly resolved by their kappa and omega values
        As the length cutoff increases, the periphery clusters decrease in size
MINMAX
    The projections are largely the same as the norm feature set projections
    At higher length cutoffs and higher perplexities, the central ball appears to separate into two subclusters
NORM
    A majority of the diverged and a significant fraction of the conserved subsequences are clustered in a large ball; however, the two classes are concentrated at opposite ends of this ball
        Though the ball is not highly structured, it is not completely uniform either
    Some generally resolved clusters of largely diverged subsequences are at the periphery of this ball
        This structure is a result of the dominance of the subsequences with an undefined kappa or omega
    With larger length cutoffs, the periphery clusters decrease in size and the ball is more clearly a continuum of conserved to diverged subsequences
ZNORM
    Similar to the norm feature set, but the boundary between the central ball and the surrounding clusters is less clear
        Within the central ball, there are more punctate subclusters as well
    Additionally, the surrounding clusters are smaller and less defined
    As in the norm feature set, the peripheral clusters decrease in size with increasing length cutoffs
These projections largely confirm the conclusion from the PCA analyses:
    Length and net charge dominate the projections of the complete feature sets, likely as a result of their scaling
    In more limited feature sets, the conserved and divergent feature sets have more overlap, but there are differences that broadly distinguish them.
        The identity of these feature is unclear
In contrast to the PCA projections, the t-SNE projections suggest subclasses of conserved and diverged subsequences
    This is a known strength of t-SNE over other methods of dimensionality reduction; it remains unclear, however, if these structures have biological relevance
Creating projections of randomly generated sequences will illuminate if these structures are inherent to the analysis or a product of selection

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""