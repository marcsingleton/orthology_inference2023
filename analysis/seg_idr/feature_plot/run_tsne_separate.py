"""Execute the tsne_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_separate.py ../sample_feats/ 3 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    The same banded fan-like structures as noted previously are observed.
        The bands are more distinct for the ordered subsequences.
MINMAX
    The same pattern as noted in the joint projections is observed except the central cluster is not composed of a gradient of the ordered and disordered subsequences.
        This reduces the number of lobes from four to two.
    The two lobes in the ordered sequences are more closely joined at the "open" end, creating a near loop. In contrast, the lobes in the disordered sequences are open and spread apart.
NORM
    The same pattern as noted in the joint projections is observed except the central cluster is not composed of a gradient of the ordered and disordered subsequences.
    Very little substructure is observed at any cutoff.
ZNORM
    The projections are nearly identical to their joint counterparts, excluding the separation of the central cluster into subclusters of the two classes.

DEPENDENCIES
../../../src/feature_plot_scripts/tsne_separate.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""