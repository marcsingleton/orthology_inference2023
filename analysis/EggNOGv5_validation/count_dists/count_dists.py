"""Plot species, duplicate, and sequence count distributions."""

import matplotlib.pyplot as plt
import os

# Compile paths of TSVs from list of directories containing TSVs and list of direct paths
dirs = ['../filter_count/7214_members/', '../filter_count/7214_noX_members/']  # Folders to check for TSVs (must end in /)
paths = ['../../../data/EggNOGv5/drosophilidae/7214_members.tsv', '../filter_unknown_realign/7214_noX_members.tsv']
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

        dist_seq = {}
        dist_species = {}
        dist_dup = {}
        for line in file:
            fields = line.rstrip().split('\t')
            species = fields[5].split(',')
            num_alignments += 1

            num_seq = int(fields[2])
            num_species = int(fields[3])
            num_dup = num_seq - num_species
            dist_seq[num_seq] = dist_seq.get(num_seq, 0) + 1
            dist_species[num_species] = dist_species.get(num_species, 0) + 1
            dist_dup[num_dup] = dist_dup.get(num_dup, 0) + 1

        # Create subdirectory to save plots and move to that directory
        if not os.path.exists(root):
            os.mkdir(root)
        os.chdir(root)

        # Print counts
        print(root.upper())
        print('number of alignments:', num_alignments)
        print()
        print('number of alignments with 10 species:', dist_species[10])
        print('fraction of alignments with 10 species:', dist_species[10] / num_alignments)
        print()
        print('number of alignments with 10 sequences:', dist_seq[10])
        print('fraction of alignments with 10 sequences:', dist_seq[10] / num_alignments)
        print()
        print('number of alignments with duplicates:', num_alignments - dist_dup[0])
        print('fraction of alignments with duplicates', (num_alignments - dist_dup[0]) / num_alignments)
        print()

        # Plot distributions
        # Distribution of number of species
        spec, spec_count = zip(*dist_species.items())
        fig1, ax1_1 = plt.subplots()
        ax1_1.bar(spec, spec_count)
        ax1_1.set_title('Distribution of Number of Species in Alignments')
        ax1_1.set_xlabel('Number of Species')
        ax1_1.set_ylabel('Count')

        ax1_2 = ax1_1.twinx()
        mn, mx = ax1_1.get_ylim()
        ax1_2.set_ylim(mn / num_alignments, mx / num_alignments)
        ax1_2.set_ylabel('Fraction')

        fig1.tight_layout()
        fig1.savefig('species_distribution.png')

        # Distribution of number of sequences
        seq, seq_count = zip(*dist_seq.items())
        fig2, ax2_1 = plt.subplots()
        ax2_1.bar(seq, seq_count, width=1, align='edge')
        ax2_1.set_title('Distribution of Number of Sequences in Alignments')
        ax2_1.set_xlabel('Number of Sequences')
        ax2_1.set_ylabel('Count')

        ax2_2 = ax2_1.twinx()
        mn, mx = ax2_1.get_ylim()
        ax2_2.set_ylim(mn / num_alignments, mx / num_alignments)
        ax2_2.set_ylabel('Fraction')

        fig2.tight_layout()
        fig2.savefig('sequence_distribution.png')

        # Distribution of number of duplicates
        seq, seq_count = zip(*dist_dup.items())
        fig3, ax3_1 = plt.subplots()
        ax3_1.bar(seq, seq_count, width=1, align='edge')
        ax3_1.set_title('Distribution of Number of Duplicates in Alignments')
        ax3_1.set_xlabel('Number of Sequences')
        ax3_1.set_ylabel('Count')

        ax3_2 = ax3_1.twinx()
        mn, mx = ax3_1.get_ylim()
        ax3_2.set_ylim(mn / num_alignments, mx / num_alignments)
        ax3_2.set_ylabel('Fraction')

        fig3.tight_layout()
        fig3.savefig('duplicate_distribution.png')

        plt.close('all')
        os.chdir('..')  # Return to initial directory

"""
OUTPUT
7214_MEMBERS
number of alignments: 13686

number of alignments with 10 species: 7551
fraction of alignments with 10 species: 0.5517316966242876

number of alignments with 10 sequences: 5352
fraction of alignments with 10 sequences: 0.3910565541429198

number of alignments with duplicates: 3700
fraction of alignments with duplicates 0.27034926201958204

7214_NOX_MEMBERS
number of alignments: 13686

number of alignments with 10 species: 7355
fraction of alignments with 10 species: 0.537410492474061

number of alignments with 10 sequences: 5298
fraction of alignments with 10 sequences: 0.38711091626479616

number of alignments with duplicates: 3628
fraction of alignments with duplicates 0.26508841151541723

7214_MEMBERS__10_10_MEMBERS
number of alignments: 5069

number of alignments with 10 species: 5069
fraction of alignments with 10 species: 1.0

number of alignments with 10 sequences: 5069
fraction of alignments with 10 sequences: 1.0

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

7214_MEMBERS__DMEL_MEMBERS
number of alignments: 11245

number of alignments with 10 species: 7551
fraction of alignments with 10 species: 0.671498443752779

number of alignments with 10 sequences: 5333
fraction of alignments with 10 sequences: 0.4742552245442419

number of alignments with duplicates: 3329
fraction of alignments with duplicates 0.29604268563806135

7214_MEMBERS__EQUAL_+5_MEMBERS
number of alignments: 7800

number of alignments with 10 species: 5069
fraction of alignments with 10 species: 0.6498717948717949

number of alignments with 10 sequences: 5069
fraction of alignments with 10 sequences: 0.6498717948717949

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

7214_NOX_MEMBERS__10_10_MEMBERS
number of alignments: 4992

number of alignments with 10 species: 4992
fraction of alignments with 10 species: 1.0

number of alignments with 10 sequences: 4992
fraction of alignments with 10 sequences: 1.0

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

7214_NOX_MEMBERS__DMEL_MEMBERS
number of alignments: 11239

number of alignments with 10 species: 7355
fraction of alignments with 10 species: 0.6544176528160869

number of alignments with 10 sequences: 5283
fraction of alignments with 10 sequences: 0.4700596138446481

number of alignments with duplicates: 3260
fraction of alignments with duplicates 0.2900613933623988

7214_NOX_MEMBERS__EQUAL_+5_MEMBERS
number of alignments: 7865

number of alignments with 10 species: 4992
fraction of alignments with 10 species: 0.6347107438016529

number of alignments with 10 sequences: 4992
fraction of alignments with 10 sequences: 0.6347107438016529

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

NOTES
Most alignments contain sequences from all ten species, but many of these likely contain duplicates.
Most alignments have no duplicates.
The distribution of alignments with duplicates appears exponential.

DEPENDENCIES
../filter_count/filter_count.py
    ../filter_count/7214_members/*.tsv
    ../filter_count/7214_noX_members/*.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/7214_noX_members.tsv
"""