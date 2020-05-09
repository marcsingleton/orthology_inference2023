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
        num_alis = 0
        ppids = set()

        # Stats relative to txids
        species_ali = {}  # Species representation in alignments
        species_seq = {}  # Species distribution in sequences

        # Stats relative to counts in alignments
        dist_nseq = {}  # Number of sequences in alignments
        dist_nspecies = {}  # Number of species in alignments
        dist_ndup = {}  # Number of duplicate species in alignments
        for line in file:
            fields = line.rstrip().split('\t')
            txids = [ppid[:4] for ppid in fields[4].split(',')]  # Gets txid
            ppids.update(fields[4].split(','))
            num_alis += 1

            # Stats relative to txids
            for txid in txids:
                species_seq[txid] = species_seq.get(txid, 0) + 1
            for txid in set(txids):
                species_ali[txid] = species_ali.get(txid, 0) + 1

            # Stats relative to counts in alignments
            nseq = len(txids)
            nspecies = len(set(txids))
            ndup = nseq - nspecies
            dist_nseq[nseq] = dist_nseq.get(nseq, 0) + 1
            dist_nspecies[nspecies] = dist_nspecies.get(nspecies, 0) + 1
            dist_ndup[ndup] = dist_ndup.get(ndup, 0) + 1

    # Create subdirectory to save plots and move to that directory
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)  # Recursive folder creation
    os.chdir('out/' + root)

    # Print counts
    print(root.upper())
    print('number of alignments:', num_alis)
    print()
    print('number of alignments with 10 species:', dist_nspecies[10])
    print('fraction of alignments with 10 species:', dist_nspecies[10] / num_alis)
    print()
    print('number of alignments with 10 sequences:', dist_nseq[10])
    print('fraction of alignments with 10 sequences:', dist_nseq[10] / num_alis)
    print()
    print('number of alignments with duplicates:', num_alis - dist_ndup[0])
    print('fraction of alignments with duplicates', (num_alis - dist_ndup[0]) / num_alis)
    print()

    # Plot stats
    # Representation of species in alignments
    labels, h_ali = zip(*sorted(species_ali.items(), key=lambda i: i[0]))
    x = list(range(1, len(labels) + 1))
    fig, ax1 = plt.subplots()
    ax1.bar(x, h_ali, tick_label=labels, align='center')
    ax1.set_xlabel('Species')
    ax1.set_ylabel('Count')
    ax1.set_title('Representation of Species in Alignments')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn / num_alis, mx / num_alis)
    ax2.set_ylabel('Fraction')

    fig.tight_layout()
    fig.savefig('species_ali_rep.png')
    plt.close()

    # Distribution of species in sequences
    labels, h_seq = zip(*sorted(species_seq.items(), key=lambda i: i[0]))
    x = list(range(1, len(labels) + 1))
    fig, ax1 = plt.subplots()
    ax1.bar(x, h_seq, tick_label=labels, align='center')
    ax1.set_xlabel('Species')
    ax1.set_ylabel('Count')
    ax1.set_title('Distribution of Species in Sequences')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn / len(ppids), mx / len(ppids))
    ax2.set_ylabel('Fraction')

    fig.tight_layout()
    fig.savefig('species_seq_dist.png')
    plt.close()

    # Correlation of alignment and sequence counts
    fig, ax = plt.subplots()
    ax.scatter(h_ali, h_seq)
    ax.set_xlabel('Species Count in Alignments')
    ax.set_ylabel('Species Count in Sequences')
    ax.set_title('Alignment and Sequence Count Correlation')

    fig.savefig('species_aliseq_corr.png')
    plt.close()

    # Distribution of number of species
    spec, spec_count = zip(*dist_nspecies.items())
    fig, ax1 = plt.subplots()
    ax1.bar(spec, spec_count)
    ax1.set_title('Distribution of Number of Species in Alignments')
    ax1.set_xlabel('Number of Species')
    ax1.set_ylabel('Count')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn / num_alis, mx / num_alis)
    ax2.set_ylabel('Fraction')

    fig.tight_layout()
    fig.savefig('num_species_dist.png')
    plt.close()

    # Distribution of number of sequences
    seq, seq_count = zip(*dist_nseq.items())
    fig, ax1 = plt.subplots()
    ax1.bar(seq, seq_count, width=1, align='edge')
    ax1.set_title('Distribution of Number of Sequences in Alignments')
    ax1.set_xlabel('Number of Sequences')
    ax1.set_ylabel('Count')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn / num_alis, mx / num_alis)
    ax2.set_ylabel('Fraction')

    fig.tight_layout()
    fig.savefig('num_sequence_dist.png')
    plt.close()

    # Distribution of number of duplicates
    seq, seq_count = zip(*dist_ndup.items())
    fig, ax1 = plt.subplots()
    ax1.bar(seq, seq_count, width=1, align='center')
    ax1.set_title('Distribution of Number of Duplicates in Alignments')
    ax1.set_xlabel('Number of Sequences')
    ax1.set_ylabel('Count')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn / num_alis, mx / num_alis)
    ax2.set_ylabel('Fraction')

    fig.tight_layout()
    fig.savefig('num_duplicate_dist.png')
    plt.close()

    os.chdir('../..')  # Return to initial directory

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

7214_MEMBERS__EQUAL_+5_MEMBERS
number of alignments: 7800

number of alignments with 10 species: 5069
fraction of alignments with 10 species: 0.6498717948717949

number of alignments with 10 sequences: 5069
fraction of alignments with 10 sequences: 0.6498717948717949

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

7214_MEMBERS__10_10_MEMBERS
number of alignments: 5069

number of alignments with 10 species: 5069
fraction of alignments with 10 species: 1.0

number of alignments with 10 sequences: 5069
fraction of alignments with 10 sequences: 1.0

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

7214_NOX_MEMBERS__EQUAL_+5_MEMBERS
number of alignments: 7865

number of alignments with 10 species: 4992
fraction of alignments with 10 species: 0.6347107438016529

number of alignments with 10 sequences: 4992
fraction of alignments with 10 sequences: 0.6347107438016529

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

7214_NOX_MEMBERS__10_10_MEMBERS
number of alignments: 4992

number of alignments with 10 species: 4992
fraction of alignments with 10 species: 1.0

number of alignments with 10 sequences: 4992
fraction of alignments with 10 sequences: 1.0

number of alignments with duplicates: 0
fraction of alignments with duplicates 0.0

DEPENDENCIES
../filter_count/filter_count.py
    ../filter_count/out/7214_members/*.tsv
    ../filter_count/out/7214_noX_members/*.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/out/7214_noX_members.tsv
"""