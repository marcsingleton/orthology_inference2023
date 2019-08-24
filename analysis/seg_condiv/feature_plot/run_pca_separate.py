"""Execute the pca_separate.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../feature_calc/ con div Conserved Diverged', shell=True)

"""
NOTES
The principal components for the conserved and diverged subsequences were analyzed separately to reveal any substructures in either that was obscured by analyzing them together
    The projections were largely the same
Nested dimensionality reductions are best for revealing substructure when classes are apparent in the primary analysis and the classes are sparse relative to the data set
    This data set fails both requirements; there are two equally-sized classes that largely overlap
    In contrast, many secondary t-SNE analyses in single cell data sets reveal sub cell types because:
        1. Clusters (cell types) were apparent in the primary reduction
        2. The number of points in each cluster is small relative to the entire data set

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../feature_calc/feature_calc.py
    ../feature_calc/features_*.tsv
"""