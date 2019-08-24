"""Execute the tsne_joint.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_plot_scripts/tsne_joint.py ../feature_calc_shuffle/ con div Conserved Diverged', shell=True)

"""
NOTES
No differences between shuffled and observed sequences for the all, len, net_charge, norm, or znorm feature sets
In the len+net_charge feature set, the shuffled sequences form clusters that are less diffuse and more punctate
    The conserved and diverged sequences still cluster, but the clusters more a mosaic of subclusters of conserved or diverged than a mixture of the sequences
    This effect is subtle, however, and largely only apparent at moderate perplexities
    

DEPENDENCIES
../../../src/feature_calc_scripts/tsne_joint.py
../feature_calc_shuffle/feature_calc_shuffle.py
    ../feature_calc_shuffle/features_?.tsv
"""