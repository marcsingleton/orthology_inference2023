"""Execute the pca_separate.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/pca_separate.py ../feature_calc/ con div Conserved Diverged', shell=True)

"""
NOTES
Differences between shuffled and observed sequences similar to those noted in run_pca_joint.py for the complete, net charge, len, and len + net charge feature sets.
The separate projections for the conserved and diverged classes highlights in the differences with regards to the presence of the kappa and omega clusters.
    The kappa and omega clusters disappear by a length cutoff of 16 in the shuffled conserved sequences.
        The clusters are still present in the original sequences until a cutoff of 32.
    In contrast, the clusters do not disppear until a length cutoff of 32 in the shuffled diverged sequences.
    In both cases, the kappa and omega clusters are denser (or more enriched) in the observed sequences.
        It is unclear if the kappa and omega clusters are more enriched in diverged than conserved sequences relative to their shuffled counterparts.
            Because the shuffled diverged sequences cluster by kappa and omega values at a higher cutoff, the amino acid pool is likely more biased towards producing sequences with no prolines or charged residues.
            It is unclear if the increased number of sequences in kappa and omega clusters outstrips their inherent propensity for the kappa and omega clusters.
            Directly counting the sequences in all categories should more clearly distinguish these possibilities.
Like in the combined projections, the separate projections differ largely in the number and the intensity of their spokes.
    The differences are diminished in the conserved sequences and exaggerated in the diverged sequences.

DEPENDENCIES
../../../src/feature_plot_scripts/pca_separate.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_*.tsv
"""