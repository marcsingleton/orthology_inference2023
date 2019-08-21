"""Execute the tsne_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_separate.py ../feature_calc/', shell=True)

"""
NOTES
Projecting the subsequences separately was successful in accentuating their substructures.
    This is consistent with the commentary on nested reductions in the PCA analysis, Clusters were apparent in the combined analysis, and these clusters were a small fraction of the overall data.
ALL
    Charge dominates the structure of the conserved projection (as a result of their unscaled length).
    Strands are also apparent in the diverged projections at the highest perplexities.
        The diverged subsequences, however, also form discrete clusters, unlike the conserved subsequences.
LEN
    Largely similar to the all feature set, as noted previously.
NET_CHARGE
    The conserved projections are largely diffuse blobs with a few smaller clusters corresponding to different values of kappa and omega.
    The diverged projections still form clusters, but they are more diffuse than those in the previous projections.
LEN+NET_CHARGE
    The clusters, particularly, those in the diverged projections are slightly more punctate.
    Elongation appears to result from feature distributions with long tails that characterize multiple clusters.
        If every cluster has a long-tailed distribution for some feature, all those clusters will attract each other which will smear their boundaries.
        The tail for the charge distribution in the conserved subsequences is so long, all the clusters are strongly attracted to a single point.
NORM
    Like in the combined projections, there are three major clusters corresponding to values of kappa and omega.
        In the conserved subsequences, these clusters are fairly uniform whereas in the diverged subsequences they display some substructure.
ZNORM
    Its star-like topology is broadly similar to the combined projection.
        However, the conserved projections are more diffuse with little to no boundaries between the central and peripherary clusters; additionally, the points are more subtle.
        In contrast, the diverged projections are more punctate along with the periphery and more clearly resolve the periphery clusters.

DEPENDENCIES
../feature_plot_scripts/tsne_separate.py
../feature_calc/feature_calc.py
    ../feature_calc/features_4.tsv
"""