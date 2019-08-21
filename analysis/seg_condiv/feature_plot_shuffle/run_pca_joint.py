"""Execute the pca_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_joint.py ../feature_calc_shuffle/', shell=True)

"""
NOTES
ALL + LEN
    Largely the same in structure.
NET_CHARGE + LEN+NET_CHARGE
    Projections are less compressed, i.e. the scale of the PC axes are smaller in the shuffled sequences.
        This is particularly apparent for the net charge feature set at length cutoffs greater than 16.
    This appears to reveal additional substructure in the shuffled sequences, but this is likely a visual artifact.
        This does suggest, however, that the random sequences have fewer outliers that increase that scale of the axes.
NORM
    At low cutoffs (< 16) the structure is nearly identical, with the kappa and omega groupings dominating the projections of both the original and shuffled sequences.
    At high cutoffs (>= 16), the kappa and omega groups largely disappear in the shuffled sequences; 
        Thus, the kappa and omega clusters are enriched in observed sequences relative to random sequences.
        This indicates that though their prominence in the projections is possibly artifactual, these clusters are not random.
    That these differences are only apparent at longer cutoffs suggests that relatively short amino acid sequences are potentially not of sufficient length to be distinguished from random sequences.
ZNORM
    No obvious differences at low cutoffs, though the shuffled sequences have fewer spokes.
    At higher cutoffs, the shuffled sequence projections are more diffuse and more clearly separated.
        These vertices or spokes in the projections thus likely correspond to non-random structures in the sequences.
        The decrease in overlap is unexplained. Are the amino acids in the diverged sequences arranged to approximate the conserved sequences (despite their different amino acid distributions)?
In conclusion, the presence of prolines or charged residues appears non-random, particularly in longer subsequences.
    Previous studies have suggested that charge (either its presence or its patterning) is a conserved feature of IDRs

DEPENDENCIES
../feature_plot_scripts/pca_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""