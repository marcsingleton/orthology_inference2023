"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_joint.py ../feature_calc_shuffle/ con div Conserved Diverged', shell=True)

"""
NOTES
ALL
    At low cutoffs, the projections are largely the same. However, the shuffled projections appear less compressed, as there are fewer outlier sequences which expand the viewing window.
    At high cutoffs, the shuffled projections maintain largely the same fan-like structure, though they appear more compressed due to the greater range of the sequences along the diagonal.
    In contrast, at high cutoffs the unshuffled projections have a parabolic structure.
MINMAX
    The observations for the norm feature set largely apply here as well up to a length cutoff of 16.
        The largest exception is that while the kappa-omega clusters are visible at lower cutoffs in the shuffled sequences, they appear distorted in substructure relative to the unshuffled sequences.
        In contrast, the unshuffled and shuffled projections have a nearly identical structure in the norm feature set at low cutoffs.
    At high cutoffs, the shuffled sequences do not separate into distinct Gaussian clusters. They are instead smeared into two distinct but overlapping bands.
        This is likely an effect of the inclusion of the length feature.
NORM
    At low cutoffs (< 16) the structure is nearly identical, with the kappa and omega groupings dominating the projections of both the original and shuffled sequences.
    At high cutoffs (>= 16), the kappa and omega groups largely disappear in the shuffled sequences; 
        Thus, the kappa and omega clusters are enriched in observed sequences relative to random sequences.
        This indicates that though their prominence in the projections is possibly artifactual, these clusters are not random.
    That these differences are only apparent at longer cutoffs suggests that relatively short amino acid sequences are potentially not of sufficient length to be distinguished from random sequences.
    Additionally, at high cutoffs (>= 16), the shuffled projections separate into two distinct but partially overlapping Gaussian clusters.
        The shuffled sequences are more distinct between classes than are the unshuffled sequences.
        Despite their different bulk amino acid compositions, within sequences they are arranged to minimize these differences.
        Likely short sequences are not of sufficient length to distinguish the classes on the basis of their amino acid biases.
ZNORM
    At low cutoffs, the shuffled sequences have fewer spokes and a kappa-omega cluster is more clearly separated.
    At higher cutoffs, the shuffled sequence projections are more diffuse and more clearly separated by conserved and diverged sequences.
        These vertices or spokes in the projections thus likely correspond to non-random structures in the sequences.
        The decrease in overlap is unexplained. Are the amino acids in the diverged sequences arranged to approximate the conserved sequences (despite their different amino acid distributions)?
In conclusion, the presence of prolines or charged residues appears non-random, particularly in longer subsequences.
    Previous studies have suggested that charge (either its presence or its patterning) is a conserved feature of IDRs

DEPENDENCIES
../../../src/feature_plot_scripts/pca_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""