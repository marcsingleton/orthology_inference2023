"""Plot boxplots of the percentage aligned for each alignment."""

import gzip
import matplotlib.pyplot as plt
import os
from Bio import AlignIO
from numpy import quantile
from numpy.random import normal


def fraction_ungapped(MSA):
    fractions = []
    length = MSA.get_alignment_length()
    for record in MSA:
        fractions.append(1 - record.seq.count('-') / length)
    return fractions

# Parameters
pathdir = [('../filter_count/7214_noX_members/equal_+5_members.tsv', '../filter_unknown_realign/align/'),
           ('../filter_unknown_realign/7214_noX_members.tsv', '../filter_unknown_realign/align/'),
           ('../../../data/EggNOGv5/drosophilidae/7214_members.tsv', '../../../data/EggNOGv5/drosophilidae/7214/')]
num_box = 15  # Number of boxes per plot
params = {'linewidths': 1, 'color': 'black', 'alpha': 0.1, 'zorder': 0}  # Boxplot display parameters

for path, dir in pathdir:
    # Create list of sorted (ID, fracs) pairs
    dists = {}
    for member in os.listdir(dir):
        with gzip.open(dir + member, 'rt') as file:
            MSA = AlignIO.read(file, 'fasta')
            ali_id = member[:5]
            dists[ali_id] = fraction_ungapped(MSA)
    sorts = sorted(dists.items(), key=lambda item: quantile(item[1], 0.5))

    # Create directory to store figures and change to that directory
    root, _ = os.path.splitext(os.path.basename(path))
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)  # Recursive folder creation
    os.chdir('out/' + root)

    # Plot analysis in groups
    for i in range(len(sorts) // num_box):
        labels, ys = zip(*sorts[num_box * i:num_box * (i + 1)])
        xs = [normal(i + 1, 0.08, len(y)) for i, y in enumerate(ys)]

        plt.figure()
        for j, (x, y) in enumerate(zip(xs, ys)):
            # Get quantiles
            q25 = quantile(y, 0.25)
            q50 = quantile(y, 0.5)
            q75 = quantile(y, 0.75)
            lw = q25 - 1.5 * (q75 - q25)
            uw = q75 + 1.5 * (q75 - q25)

            # Plot lines and markers
            plt.hlines(q25, j + 0.75, j + 1.25, **params)
            plt.hlines(q75, j + 0.75, j + 1.25, **params)
            plt.hlines(lw, j + 0.9375, j + 1.0625, **params)
            plt.hlines(uw, j + 0.9375, j + 1.0625, **params)
            plt.vlines(j + 0.75, q25, q75, **params)
            plt.vlines(j + 1.25, q25, q75, **params)
            plt.vlines(j + 1, lw, q25, **params)
            plt.vlines(j + 1, q75, uw, **params)
            plt.hlines(q50, j + 0.75, j + 1.25, color='orange', zorder=0)
            plt.scatter(x, y, s=10, color='black', edgecolors='none', alpha='0.5')

        # Format axes
        plt.xlabel('Alignment ID')
        plt.xticks(list(range(1, j + 2)), labels, rotation=60)
        plt.ylabel('Ungapped Percentage')
        plt.ylim((-0.1, 1.1))
        plt.subplots_adjust(bottom=0.2)

        plt.savefig(f'boxplot{i}.png')
        plt.close()

    os.chdir('../..')

"""
NOTES
Using boxplots to identify "outlier" alignments is likely not a good strategy. For example, some sequences are so highly
conserved that the IQR is nearly 0, which excludes sequences of reasonable variability. In other cases, the distribution
of alignment percentages form clusters (likely the result of actual clusters in sequence space). If the clusters are on
either side of the median, the IQR is inflated, which can result in the inclusion of clearly spurious sequences. In
general this approach is flawed because it does not account for relationships between sequences, and attempts to remove
spurious sequences by a bulk measure. A phylogenetic probabilitistic approach would be a more precise method of
identifying highly divergent sequences.

DEPENDENCIES
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/7214_noX_members.tsv
    ../filter_unknown_realign/align/
../filter_count/filter_count.py
    ../filter_count/7214_noX_members/equal_+5_members.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
    ../../../data/EggNOGv5/drosophilidae/7214/
"""