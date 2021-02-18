"""Trim alignments, removing spurious sequences and regions."""

import os
import re
from itertools import chain
from math import ceil, exp, log

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import skbio
from matplotlib.gridspec import GridSpec
from src.draw import draw_msa


def plot_msa_lines(msa, lines, figsize=(15, 7.5),
                   msa_height=1, data_height=1, hspace=0.75, sym_length=7, sym_height=7,
                   lines_min=None, lines_max=None,
                   block_cols=None, aa2color=None):
    # Define functions and globals
    ROWS, COLS = len(msa), len(msa[0])
    RATIO = figsize[0] / figsize[1]

    def get_dims(block_cols):
        plot_length = block_cols  # Length of final plot
        block_num = COLS // block_cols - (1 if COLS % block_cols == 0 else 0)  # Number of blocks in addition to the first
        plot_height = 2 * (ROWS + hspace) * block_num + (2 + hspace) * ROWS  # Height of final image
        return plot_length, plot_height

    def get_aspect(block_cols):
        plot_length, plot_height = get_dims(block_cols)
        return plot_length / plot_height

    def get_block_cols():
        # Use binary search to find interval containing optimal block_cols
        interval = (1, COLS)
        while interval[1] - interval[0] > 1:
            i1 = (interval[0], (interval[0] + interval[1]) // 2)
            i2 = ((interval[0] + interval[1]) // 2, interval[1])
            if (get_aspect(i1[0]) - RATIO) * (get_aspect(i1[1]) - RATIO) < 0:
                interval = i1
            elif (get_aspect(i2[0]) - RATIO) * (get_aspect(i2[1]) - RATIO) < 0:
                interval = i2
            else:
                break
        block_cols = min(interval, key=lambda x: abs(get_aspect(x) - RATIO))  # Choose value that minimizes difference

        # Ensure last block is at least 50% of block_cols
        if COLS % block_cols < 0.5 * block_cols:  # Guarantees at least two blocks
            blocks_im = COLS // block_cols  # Total blocks minus 1
            block_cols += ceil((COLS % block_cols) / blocks_im)  # Distribute excess to other blocks
        return block_cols

    # Set options
    if block_cols is None:
        block_cols = get_block_cols()
        block_num = COLS // block_cols + (1 if COLS % block_cols > 0 else 0)  # Number of blocks
    if isinstance(lines, list):
        lines = np.array(lines)
    if lines.ndim == 1:
        lines = np.expand_dims(lines, axis=0)
    if lines_min is None:
        lines_min = lines.min() - 0.05 * (lines.max() - lines.min())
    if lines_max is None:
        lines_max = lines.max() + 0.05 * (lines.max() - lines.min())
    if lines_min == lines_max:
        lines_min -= 0.5
        lines_max += 0.5
    block_rows = len(msa)

    im = draw_msa(msa, im_cols=len(msa[0]),
                  sym_length=sym_length, sym_height=sym_height, aa2color=aa2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2*block_num, 1, figure=fig,
                  height_ratios=[msa_height if i % 2 == 0 else data_height for i in range(2*block_num)],
                  hspace=hspace)
    for i in range(block_num):
        msa_ax = fig.add_subplot(gs[2*i:2*i+1, :])
        lines_ax = fig.add_subplot(gs[2*i+1:2*(i+1), :], sharex=msa_ax, aspect=block_rows*data_height/(msa_height * (lines_max - lines_min)))

        block = im[:, i*sym_length*block_cols:(i+1)*sym_length*block_cols]
        msa_ax.imshow(block, extent=[i*block_cols, i*block_cols + block.shape[1]//sym_length, 0, block_rows], origin='lower')
        msa_ax.axis('off')

        for line in lines:
            lines_ax.plot(list(range(i*block_cols, i*block_cols + block.shape[1]//sym_length)),
                          line[i * block_cols:i * block_cols + block.shape[1] // sym_length])
        lines_ax.set_ylim(lines_min, lines_max)
        lines_ax.tick_params(labelsize=8)
    return fig


def plot_msa(msa, figsize=(15, 7.5),
            labels=None, labelsize=6,
            hspace=0.5, sym_length=7, sym_height=7,
            block_cols=None, aa2color=None):
    # Define functions and globals
    ROWS, COLS = len(msa), len(msa[0])
    RATIO = figsize[0] / figsize[1]

    def get_dims(block_cols):
        plot_length = block_cols  # Length of final plot
        block_num = COLS // block_cols - (1 if COLS % block_cols == 0 else 0)  # Number of blocks in addition to the first
        plot_height = (ROWS + hspace) * block_num + (1 + hspace) * ROWS  # Height of final image
        return plot_length, plot_height

    def get_aspect(block_cols):
        plot_length, plot_height = get_dims(block_cols)
        return plot_length / plot_height

    def get_block_cols():
        # Use binary search to find interval containing optimal block_cols
        interval = (1, COLS)
        while interval[1] - interval[0] > 1:
            i1 = (interval[0], (interval[0] + interval[1]) // 2)
            i2 = ((interval[0] + interval[1]) // 2, interval[1])
            if (get_aspect(i1[0]) - RATIO) * (get_aspect(i1[1]) - RATIO) < 0:
                interval = i1
            elif (get_aspect(i2[0]) - RATIO) * (get_aspect(i2[1]) - RATIO) < 0:
                interval = i2
            else:
                break
        block_cols = min(interval, key=lambda x: abs(get_aspect(x) - RATIO))  # Choose value that minimizes difference

        # Ensure last block is at least 50% of block_cols
        if COLS % block_cols < 0.5 * block_cols:  # Guarantees at least two blocks
            blocks_im = COLS // block_cols  # Total blocks minus 1
            block_cols += ceil((COLS % block_cols) / blocks_im)  # Distribute excess to other blocks
        return block_cols

    # Set options
    if labels is None:
        labels = []
    if block_cols is None:
        block_cols = get_block_cols()
        block_num = COLS // block_cols + (1 if COLS % block_cols > 0 else 0)  # Number of blocks
    block_rows = len(msa)

    im = draw_msa(msa, im_cols=len(msa[0]),
                  sym_length=sym_length, sym_height=sym_height, aa2color=aa2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(block_num, 1, figure=fig,
                  hspace=hspace)
    for i in range(block_num):
        msa_ax = fig.add_subplot(gs[i:i+1, :])

        block = im[:, i*sym_length*block_cols:(i+1)*sym_length*block_cols]
        msa_ax.imshow(block, extent=[i*block_cols, i*block_cols + block.shape[1]//sym_length, 0, block_rows], origin='lower')
        msa_ax.set_yticks([x+0.5 for x in range(ROWS)])
        msa_ax.set_yticklabels(labels)
        msa_ax.tick_params(axis='y', length=0, labelsize=labelsize)
        msa_ax.spines['bottom'].set_visible(False)
        msa_ax.spines['left'].set_visible(False)
    return fig


def get_segments(msa, region):
    block = msa.loc[:, region]
    segments = {i: {'region': region, 'index': i, 'slices': [], 'support': 0} for i in range(len(block))}
    starts = {i: None for i in range(len(block))}  # Slice starts
    for j, col in enumerate(block.iter_positions()):
        support = len(col) - str(col).count('-') - 1
        for i, sym in enumerate(str(col)):
            if sym != '-':
                segments[i]['support'] += support
                if starts[i] is None:
                    starts[i] = j
            elif starts[i] is not None:
                segments[i]['slices'].append(slice(starts[i] + region.start, j + region.start))
                starts[i] = None
    for i, start in starts.items():
        if start is not None:
            segments[i]['slices'].append(slice(starts[i] + region.start, region.stop))
    return [segment for segment in segments.values() if segment['slices']]


def is_trimmed(slices, support, gap_score, gap_bias):
    b = 0
    w1, w2, w3, w4 = 1, -1, -1, 2
    threshold = 0.5

    length = sum([s.stop - s.start for s in slices])
    e = b + w1*length ** 2 + w2*support + w3*gap_score + w4*gap_bias
    p = 1 / (1 + exp(-e))
    return p > threshold


# TAXONOMIC CONDITIONS
conditions = [(set(['dinn', 'dgri']), 1),
              (set(['dnov', 'dvir', 'dhyd', 'dmoj', 'dnav']), 4),
              (set(['dper', 'dpse']), 1),
              (set(['dobs', 'dgua', 'dsub']), 2),
              (set(['dana', 'dbip', 'dkik', 'dser']), 3),
              (set(['dele', 'dfic', 'dtak', 'dbia']), 3),
              (set(['deug', 'dere']), 1),
              (set(['dmel']), 1),
              (set(['dmau', 'dsec']), 1)]

# CONSERVED REGIONS PARAMETERS
c_fraction = 0.15  # Maximum gap fraction in conserved columns
c_close = 3  # Size of closing element
c_length = 15  # Minimum number of columns in conserved regions

# CONSERVED REGIONS TRIMMING PARAMETERS
tc_lambda = 0.15  # Decay rate of trim signal
tc_threshold = 5  # Minimum number of residues
tc_cutoff = tc_threshold / 1000  # Minimum increment at which to stop propagating signal

# GAP REGIONS PARAMETERS
g_fraction = 0.85  # Minimum gap fraction in gap columns
l_sigma = 2  # Filter size for calculation of local gap bias
nl_lambda = 0.05  # Decay rate of nonlocal gap bias
nl_cutoff = 0.001  # Minimum increment at which to stop propagating signal

# GAP REGIONS TRIMMING PARAMETERS
tg_lambda = 0.75  # Decay rate of trim signal
tg_threshold = 1  # Minimum number of residues
tg_cutoff = tg_threshold / 1000  # Minimum increment at which to stop propagating signal

if not os.path.exists('out/'):
    os.mkdir('out/')

paths1 = [f'../align_fastas1/out/{file}' for file in os.listdir('../align_fastas1/out/') if file.endswith('.mfa')]
paths2 = [f'../align_fastas2-2/out/{file}' for file in os.listdir('../align_fastas2-2/out/') if file.endswith('.mfa')]
for path in chain(paths1, paths2):
    # 1 Load data and perform checks
    # 1.1 Load MSA
    msa1 = skbio.read(path, format='fasta', into=skbio.TabularMSA, constructor=skbio.Protein)
    for seq in msa1:
        match = re.match('ppid=([A-Za-z_0-9]+)\|gnid=([A-Za-z_0-9]+)\|spid=([a-z]+)', seq.metadata['id'])
        d = {key: value for key, value in zip(['ppid', 'gnid', 'spid'], match.groups())}
        seq.metadata.update(d)

    # 1.2 Test for species groups
    spids = set([seq.metadata['spid'] for seq in msa1])
    if any([len(spids & group) < num for group, num in conditions]):
        continue

    # 1.3 Create list of symbols representing final trimmed alignment
    syms_list = [list(str(seq)) for seq in msa1]

    # 2.1 Calculate gap scores
    scores = np.zeros(msa1.shape[1])
    for i, col in enumerate(msa1.iter_positions()):
        scores[i] = col.count('-')

    # 2.2 Get conserved regions
    mask = ndimage.label(ndimage.binary_closing(scores / len(msa1) < c_fraction), structure=c_close * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= c_length]
    msa2 = msa1.loc[:, regions]
    trim2full, i = {}, 0  # Trimmed to full MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            trim2full[i] = j
            i += 1

    # 2.3 Trim segments near gaps
    for i, seq in enumerate(msa2):
        # 2.3.1 Get gap segments
        gaps1 = np.array([sym == '-' for sym in str(seq)])  # Sequence coded as boolean gap or non-gap
        trim_signal = np.zeros(len(gaps1))  # Deletion signals
        for region, in ndimage.find_objects(ndimage.label(gaps1)[0]):
            length = region.stop - region.start
            max_k = -log(tc_cutoff / length) / tc_lambda

            k = 0
            while region.start - k - 1 >= 0 and k < max_k:
                if not gaps1[region.start - k - 1]:
                    v = length * exp(-tc_lambda*k)
                    trim_signal[region.start - k - 1] += v
                k += 1

            k = 0
            while region.stop + k <= len(trim_signal) - 1 and k < max_k:
                if not gaps1[region.stop + k]:
                    v = length * exp(-tc_lambda*k)
                    trim_signal[region.stop + k] += v
                k += 1

        # 2.3.2 Trim non-gap segments
        gaps2 = trim_signal > tc_threshold  # Sequence coded as boolean gap or non-gap after signal propagation
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            for j in range(region.start, region.stop):  # Iterate over positions to not delete segments between boundaries
                syms[trim2full[j]] = '-'

    # 3.1 Get gap regions
    mask = ndimage.label(scores / len(msa1) > g_fraction)[0]
    regions = [region for region, in ndimage.find_objects(mask)]

    # 3.2 Get segments in regions
    segments = []
    for region in regions:
        segments.extend(get_segments(msa1, region))
    segments = sorted(segments, key=lambda x: sum([s.stop-s.start for s in x['slices']]), reverse=True)

    # 3.3 Get local gap bias
    # 3.3.1 Make MSA of non-gap regions
    mask = ndimage.label(mask == 0)[0]  # Invert previous mask
    regions = [region for region, in ndimage.find_objects(mask)]
    msa3 = msa1.loc[:, regions]
    full2trim, i = {}, 0  # Full to trimmed MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            full2trim[j] = i
            i += 1

    # 3.3.2 Calculate gap scores of non-gap MSA
    scores = np.zeros(msa3.shape[1])
    for i, col in enumerate(msa3.iter_positions()):
        scores[i] = col.count('-')
    local_signal = ndimage.gaussian_filter1d(scores / len(msa3), sigma=l_sigma, mode='constant', cval=1)

    # 3.4 Make signals arrays (for nonlocal gap bias and trimming near gaps)
    gap_signals = np.zeros(msa1.shape)
    trim_signals = np.zeros(msa1.shape)
    gaps_array1 = np.array([[sym == '-' for sym in str(seq)] for seq in msa1])  # Sequences coded as boolean gap or non-gap

    # 3.5 Get trim slices
    trims1 = []
    for segment in segments:
        # 3.5.1 Unpack names
        region = segment['region']
        index = segment['index']
        slices = segment['slices']
        length = sum([s.stop-s.start for s in slices])

        # 3.5.2 Get local gap bias
        if region.start == 0:
            signal1 = 1
        else:
            idx1 = full2trim[region.start-1]
            signal1 = local_signal[idx1]
        if region.stop == msa1.shape[1]:
            signal2 = 1
        else:
            idx2 = full2trim[region.stop]
            signal2 = local_signal[idx2]
        local_bias = (signal1 + signal2) / 2

        # 3.5.3 Get nonlocal gap bias
        gap_signal = gap_signals[index]
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        nonlocal_bias = gap_signal[start:stop].sum()

        # 3.5.4 Update trims, nonlocal bias, and deletion signal
        if is_trimmed(slices, segment['support'], local_bias, nonlocal_bias):
            trims1.extend(slices)

            # Nonlocal bias
            max_k = -log(nl_cutoff / length) / nl_lambda

            k = 0
            while start - k - 1 >= 0 and k < max_k:
                if not gaps_array1[index, start - k - 1]:
                    v = length * exp(-nl_lambda*k)
                    gap_signal[start - k - 1] += v
                k += 1

            k = 0
            while stop + k <= len(gap_signal) - 1 and k < max_k:
                if not gaps_array1[index, stop + k]:
                    v = length * exp(-nl_lambda*k)
                    gap_signal[stop + k] += v
                k += 1

            # Trim signal
            max_k = -log(tg_cutoff / length) / tg_lambda

            k = 0
            while start - k - 1 >= 0 and k < max_k:
                if not gaps_array1[index, start - k - 1]:
                    v = length * exp(-tg_lambda*k)
                    trim_signals[index, start - k - 1] += v
                k += 1

            k = 0
            while stop + k <= len(trim_signals[index]) - 1 and k < max_k:
                if not gaps_array1[index, stop + k]:
                    v = length * exp(-tg_lambda*k)
                    trim_signals[index, stop + k] += v
                k += 1
    trims1 = sorted(trims1, key=lambda x: (x.start, x.start-x.stop))  # Start then -length

    # 3.6 Merge and invert slices
    # 3.6.1 Merge slices
    trims2 = []
    if trims1:
        start, stop = trims1[0].start, trims1[0].stop
        for trim in trims1:
            if trim.start > stop:
                trims2.append(slice(start, stop))
                start = trim.start
            if trim.stop > stop:
                stop = trim.stop
        trims2.append(slice(start, stop))

    # 3.6.2 Invert slices
    trims3 = []
    start = 0
    for trim in trims2:
        trims3.append(slice(start, trim.start))
        start = trim.stop
    trims3.append(slice(start, msa3.shape[1]))

    # 3.7 Combine trims (segments and columns) to yield final alignment
    gaps_array2 = trim_signals > tg_threshold
    for i, gaps2 in enumerate(gaps_array2):
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            syms[region.start:region.stop] = (region.stop - region.start) * ['-']
    seqs = [skbio.Protein(''.join(syms), metadata=seq.metadata) for syms, seq in zip(syms_list, msa1)]
    msa4 = skbio.TabularMSA(seqs).loc[:, trims3]
    skbio.write(msa4, format='fasta', into=f'out/{path[-8:]}', max_width=80)

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
"""