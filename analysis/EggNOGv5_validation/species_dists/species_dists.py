"""Plot distribution of species in alignments and all sequences."""

import matplotlib.pyplot as plt
import os

# Compile paths of TSVs from list of directories containing TSVs and list of direct paths
dirs = ['../filter_count/out/7214_members/', '../filter_count/out/7214_noX_members/']  # Folders to check for TSVs (must end in /)
paths = ['../../../data/EggNOGv5/drosophilidae/7214_members.tsv', '../filter_unknown_realign/out/7214_noX_members.tsv']
for dir in dirs:
    paths.extend([dir + file for file in os.listdir(dir) if file.endswith('.tsv')])

# Create names of output directories and rename collisions
roots = [os.path.splitext(os.path.basename(path))[0] for path in paths]
roots_idx = {}
for i, root in enumerate(roots):
    try:
        roots_idx[root].append(i)
    except KeyError:
        roots_idx[root] = [i]
for root, idx in roots_idx.items():
    if len(idx) > 1:
        for i in idx:
            pieces = os.path.abspath(paths[i]).split(os.sep)
            roots[i] = pieces[-2] + '__' + root

for path, root in zip(paths, roots):
    with open(path) as file:
        num_alignments = 0
        num_seq = 0
        dist_alignments = {}
        dist_seq = {}
        for line in file:
            fields = line.rstrip().split('\t')
            species = [protein[0:4] for protein in fields[4].split(',')]  # Slice first 4 characters of proteins in alignment
            num_alignments += 1
            for spec in species:
                num_seq += 1
                dist_seq[spec] = dist_seq.get(spec, 0) + 1
            for spec in set(species):
                dist_alignments[spec] = dist_alignments.get(spec, 0) + 1

    # Create subdirectory to save plots and move to that directory
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)  # Recursive folder creation
    os.chdir('out/' + root)

    # Plot distributions
    # Representation of species in alignments
    labels1, h1 = zip(*sorted(dist_alignments.items(), key=lambda i: i[0]))
    x1 = list(range(1, len(labels1) + 1))
    fig1, ax1_1 = plt.subplots()
    ax1_1.bar(x1, h1, tick_label=labels1, align='center')
    ax1_1.set_xlabel('Species')
    ax1_1.set_ylabel('Count')
    ax1_1.set_title('Representation of Species in Alignments')

    ax1_2 = ax1_1.twinx()
    mn, mx = ax1_1.get_ylim()
    ax1_2.set_ylim(mn / num_alignments, mx / num_alignments)
    ax1_2.set_ylabel('Fraction')

    fig1.tight_layout()
    fig1.savefig('alignment_representation.png')

    # Distribution of species in sequences
    labels2, h2 = zip(*sorted(dist_seq.items(), key=lambda i: i[0]))
    x2 = list(range(1, len(labels2) + 1))
    fig2, ax2_1 = plt.subplots()
    ax2_1.bar(x2, h2, tick_label=labels2, align='center')
    ax2_1.set_xlabel('Species')
    ax2_1.set_ylabel('Count')
    ax2_1.set_title('Distribution of Species in Sequences')

    ax2_2 = ax2_1.twinx()
    mn, mx = ax2_1.get_ylim()
    ax2_2.set_ylim(mn / num_seq, mx / num_seq)
    ax2_2.set_ylabel('Fraction')

    fig2.tight_layout()
    fig2.savefig('sequence_distribution.png')

    # Correlation of alignment and sequence counts
    fig3, ax3 = plt.subplots()
    ax3.scatter(h1, h2)
    ax3.set_xlabel('Species Count in Alignments')
    ax3.set_ylabel('Species Count in Sequences')
    ax3.set_title('Alignment and Sequence Count Correlation')

    fig3.savefig('correlation.png')

    plt.close('all')
    os.chdir('../..')  # Return to initial directory

"""
DEPENDENCIES
../filter_count/filter_count.py
    ../filter_count/out/7214_members/*.tsv
    ../filter_count/out/7214_noX_members/*.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/out/7214_noX_members.tsv
"""