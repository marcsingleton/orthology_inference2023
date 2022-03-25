"""Calculate nucleotide frequencies in alignments."""

import os

import matplotlib.pyplot as plt
import skbio


def grouped_bar(groups, group_width, bar_width, file_label, bar_labels=None, bar_colors=None):
    group_labels = sorted(groups)
    if bar_labels is None:
        bar_labels = list({bar_label for group in groups.values() for bar_label in group})
    if bar_colors is None:
        bar_colors = [f'C{i%10}' for i in range(len(bar_labels))]

    xs = [i * group_width for i in range(len(groups))]
    lim = len(bar_labels) // 2
    shift = 0 if len(bar_labels) % 2 == 1 else 0.5
    dxs = [bar_width * (x + shift) for x in range(-lim, len(bar_labels) - lim)]

    plt.figure(figsize=(8, 4))
    for dx, bar_label, bar_color in zip(dxs, bar_labels, bar_colors):
        plt.bar([x + dx for x in xs], [groups[group_label].get(bar_label, 0) for group_label in group_labels],
                width=bar_width, color=bar_color, label=bar_label)
    plt.xticks(xs, group_labels, rotation=60, fontsize='small')
    plt.xlabel('Species')
    plt.ylabel('Frequency')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 'small'})
    plt.subplots_adjust(bottom=0.2, left=0.1, right=0.9)
    plt.savefig(f'out/{file_label}.png')
    plt.close()


# Count nucleotides
counts = {}
for path in [path for path in os.listdir('../align_AA2NT/out/') if path.endswith('.afa')]:
    msa = skbio.read(f'../align_AA2NT/out/{path}', 'fasta')
    for seq in msa:
        spid = seq.metadata['id'][-4:]
        for sym in str(seq):

            if sym != '-':
                try:
                    counts[spid][sym] = counts[spid].get(sym, 0) + 1
                except KeyError:
                    counts[spid] = {sym: 1}

if not os.path.exists('out/'):
    os.mkdir('out/')

# Plot nucleotide frequencies
freqs = {}
for spid, count in counts.items():
    total = sum(count.values())
    freqs[spid] = {sym: num / total for sym, num in count.items()}

grouped_bar(freqs, 1.5, 0.25, 'nt_freqs', ['A', 'T', 'G', 'C'], ['C0', 'C1', 'C3', 'C2'])

# Plot AT/GC frequencies
ATGCs = {}
for spid, freq in freqs.items():
    ATGCs[spid] = {'AT': freq['A'] + freq['T'], 'GC': freq['G'] + freq['C']}

grouped_bar(ATGCs, 1, 0.3, 'ATGC_freqs', ['AT', 'GC'])

"""
DEPENDENCIES
../align_AA2NT/align_AA2NT.py
    ../align_AA2NT/out/*.afa
"""