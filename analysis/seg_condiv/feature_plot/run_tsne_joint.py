"""Execute the tsne_joint.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc/', shell=True)

"""
NOTES
The perplexity has a significant impact on the topology of the projection.
    At the lowest perplexity, the projections are always "balls" of uniform density; though, conserved and diverged subsequences form distinct regions.
    Significant structures do not tend to emerge until the perplexity is 20.
    The projections at the highest perplexity are usually highly structured. None appear to have a "ball" topology, indicating the perplexity is not set too high.
ALL
    The diverged subsequences form elongated "threads"; these are likely clustered by length.
    The conserved subsequences form a series of clusters along the periphery.
        Many are small, on the order of tens of points, and poorly resolved from adjacent clusters.
        Some are larger, but these often have substructure as well suggesting even broadly similar sequences are composed of heterogeneous subpopulations.
        Though the size and poor resolution of the clusters are poor, they are largely more defined than simulations with multidimensional Gaussians.
    The conserved and diverged subsequences are not completely resolved. In particular, though they are largely composed of conserved subsequences, the periphery clusters also contain some diverged subsequences.
    The periphery clusters are highly resolved by their kappa and omega values.
LEN
    The projections are largely the same as those from the complete feature set.
NET_CHARGE
    These projections do not contain the elongated "threads", indicating they were a result of net charge rather than length.
        Now the diverged subsequences are a ball in the center with little to no substructure; though at higher perplexities, some structure without distinct clusters is visible.
    The clusters of the conserved subsequences along the periphery are largely the same as in other feature sets.
        The central and periphery clusters are also highly resolved by kappa and omega values.
LEN+NET_CHARGE
    Similar to the feature set without net charge, but the conserved and diverged subsequences as well as individual clusters are less resolved at higher perplexities.
NORM
    A majority of the diverged and a significant fraction of the conserved subsequences are clustered in a large ball; however, the two classes are concentrated at opposite ends of this ball.
    Some generally resolved clusters of largely diverged subsequences are at the peripherary of this ball.
        This structure is a result of the dominance of the subsequences with an undefined kappa or omega.
ZNORM
    Similar to the norm feature set, but the boundary between the central ball and the surrounding clusters is less clear.
    Additionally, the surrounding clusters are smaller and less defined.
These projections largely confirm the conclusion from the PCA analyses:
    Length and net charge dominate the projections of the complete feature sets, likely as a result of their scaling.
    In more limited feature sets, the conserved and divergent feature sets have more overlap, but there are differences that broadly distinguish them.
        The identity of these feature is unclear.
In contrast to the PCA projections, the t-SNE projections suggest subclasses of conserved and diverged subsequences.
    This is a known strength of t-SNE over other methods of dimensionality reduction; it remains unclear, however, if these structures have biological relevance.
Creating projections of randomly generated sequences will illuminate if these structures are inherent to the analysis or a product of selection.

DEPENDENCIES
../feature_calc_scripts/tsne_joint.py
../feature_calc/feature_calc.py
    ../feature_calc/features_4.tsv
"""