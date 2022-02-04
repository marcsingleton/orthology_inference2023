"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/pca_joint.py ../feature_calc_shuffle/out/ 1 ordered Ordered Disordered', shell=True)

"""
NOTES
ALL
    The unshuffled projections are parabolic whereas the shuffled projections are fan-shaped, as observed in other PCAs.
        This is likely an effect of the shuffling significantly reducing the spread of the SCD.
        Examinations of the model summaries show a decrease in the contribution of the SCD to the overall variance.
        Thus, we are likely viewing a projection along a different axis.
MINMAX
    At low cutoffs. the projections are nearly identical, with distorted kappa-omega clusters.
    At high cutoffs (>= 16), the projections differ significantly.
        The kappa-omega clusters disappear for both, and the central cluster is "smeared" relative to the corresponding cluster in the norm feature set.
        However, in the observed sequences the classes largely overlap (excluding a concentration of ordered subsequences in the center) whereas in the shuffled sequences they form a clear gradient.
            This is likely an effect of the decrease in the variance of certain components changing the "viewing perspective" as discussed in the notes for the norm feature set.
                For example, PC2 in the shuffled subsequences contains a significant contribution from loglen whereas neither PC1 nor PC2 contain loglen in their first five components.
                Thus, this structure likely exists in the observed sequences since the loglen distribution is identical by construction, but we only observe it in the shuffled sequences since its relative contribution to the variance is much greater.
NORM
    At low cutoffs, the kappa-omega clusters are visible in both projections.
        However, the shuffled sequences show more separation between the ordered and disordered sequences within the central cluster.
    The kappa-omega clusters persist in the observed sequences until a cutoff of 16 whereas they disappear in the shuffled sequences by a cutoff of 8.
    As the cutoff increases, the separation between the ordered and disordered subsequences increases, but the separation is greater in the shuffled subsequences.
        Certainly for both subsequences, longer cutoffs allow the amino acid biases present in both to increase in importance since with more information a subsequence can be assigned to a category with higher confidence.
        However, it appears that these differences are a greater proportion of the overall variance when the sequences are shuffled.
        Potentially this structure is already present in the unshuffled sequences, but because other features contribute more to the variance due to the non-random distribution of amino acids in sequences, we are viewing the data along an axis that is less biologically-informative.
        For example, the charge distributions often have outliers due to runs of charges, but this can occur in both ordered and disordered subsequences, so this axis is not necessarily the most informative for finding differences even though it is an axis that significantly contributes to the total variance.    
ZNORM
    Even at the lowest cutoffs, there is separation between the ordered and disordered sequences, and by a cutoff of 8, this differences is significant.
    The observed projections have more flattened and angular projections whereas the shuffled projections are more symmetric and Gaussian in structure.
The differences observed in the projections are potentially not the result of large differences in their structures, but more an artifact of the projections.
    This effect was noted in the dramatic changes in the appearance of projections in the drop feature sets.
In future experiments, the data should be displayed using "cross" transformation matrices, i.e. projecting the shuffled data with the observed PCA and vice versa.
    This will normalize for the effect of the analyses having dramatically different "views" of the data.

DEPENDENCIES
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/out/features_*.tsv
"""