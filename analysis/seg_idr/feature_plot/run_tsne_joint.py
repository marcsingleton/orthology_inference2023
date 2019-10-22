"""Execute the tsne_joint.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../sample_feats/ 3 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    At all length cutoffs, the projections form banded fan-like structures.
    At short length cutoffs, periphery clusters are visible. These correspond to kappa-omega clusters.
MINMAX
    Like the norm feature set, except the central cluster forms two lobes.
        At short cutoffs these lobes are joined whereas at long cutoffs they separate.
        At long cutoffs the gradient across the two classes for each lobe sharpens, forming two lobes within each major lobe.
            Thus, the final structure is four-lobed
    The central cluster retains significant substructure at moderate cutoffs, though it is significantly reduced at the highest levels.
NORM
    At short length cutoffs, there is one large central cluster with two periphery kappa-omega clusters.
        The ordered and disordered classes form a gradient across this central cluster.
    At long length cutoffs, the kappa-omega clusters disappear and the gradient sharpens, forming two distinct clusters at the highest cutoff.
ZNORM
    At short length cutoffs, the projections have the characteristic "vertices" with periphery kappa-omega clusters.
        The ordered subsequences form a nucleus at the center of this cluster surrounded by disordered subsequences.
    At long length cutoffs, the kappa-omega clusters disappear and the vertices are significantly reduced in number and sharpness (but they are not absent).
        The ordered and disordered subsequences form nearly separate subclusters at opposite sides of the central cluster.

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../sample_feats/run_sample_feats.py
    ../sample_feats/features_*.tsv
"""