"""Trim columns and segments from an alignment."""

from math import exp, log

import numpy as np
import scipy.ndimage as ndimage
import skbio

# CONSERVED REGIONS PARAMETERS
c_fraction = 0.15  # Maximum gap fraction in conserved columns
c_close = 3  # Size of closing element
c_length = 15  # Minimum number of columns in conserved regions

# CONSERVED REGIONS TRIMMING PARAMETERS
tc_rate = 0.15  # Decay rate of trim signal
tc_threshold = 5  # Minimum number of residues
tc_cutoff = tc_threshold / 1000  # Minimum increment at which to stop propagating signal

# GAP REGIONS PARAMETERS
g_fraction = 0.85  # Minimum gap fraction in gap columns
l_sigma = 2  # Filter size for calculation of local gap bias
nl_rate = 0.05  # Decay rate of nonlocal gap bias
nl_cutoff = 0.001  # Minimum increment at which to stop propagating signal

# GAP REGIONS TRIMMING PARAMETERS
tg_rate = 0.75  # Decay rate of trim signal
tg_threshold = 1  # Minimum number of residues
tg_cutoff = tg_threshold / 1000  # Minimum increment at which to stop propagating signal


def trim_msa(msa):
    """Trim MSA by removing large insertions and segments near indels."""
    msa1 = msa  # Re-name to accord with following naming convention
    syms_list = [list(str(seq)) for seq in msa1]   # Create list of symbols representing final trimmed alignment

    # 1.1 Calculate gap scores
    scores = np.zeros(msa1.shape[1])
    for i, col in enumerate(msa1.iter_positions()):
        scores[i] = col.count('-')

    # 1.2 Get conserved regions
    mask = ndimage.label(ndimage.binary_closing(scores / len(msa1) < c_fraction), structure=c_close * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= c_length]
    msa2 = msa1.loc[:, regions]
    trim2full, i = {}, 0  # Trimmed to full MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            trim2full[i] = j
            i += 1

    # 1.3 Trim segments near gaps
    for i, seq in enumerate(msa2):
        # 1.3.1 Get gap segments
        gaps1 = np.array([sym == '-' for sym in str(seq)])  # Sequence coded as boolean gap or non-gap
        trim_signal = np.zeros(len(gaps1))  # Deletion signals
        for region, in ndimage.find_objects(ndimage.label(gaps1)[0]):
            length = region.stop - region.start
            propagate(region.start, region.stop, length, trim_signal, gaps1, tc_rate, tc_cutoff)

        # 1.3.2 Trim non-gap segments
        gaps2 = trim_signal > tc_threshold  # Sequence coded as boolean gap or non-gap after signal propagation
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            for j in range(region.start, region.stop):  # Iterate over positions to not delete segments between boundaries
                syms[trim2full[j]] = '-'

    # 2.1 Get gap regions
    mask = ndimage.label(scores / len(msa1) > g_fraction)[0]
    regions = [region for region, in ndimage.find_objects(mask)]

    # 2.2 Get segments in regions
    segments = []
    for region in regions:
        segments.extend(get_segments(msa1, region))
    segments = sorted(segments, key=lambda x: sum([s.stop-s.start for s in x['slices']]), reverse=True)

    # 2.3 Get local gap bias
    # 2.3.1 Make MSA of non-gap regions
    mask = ndimage.label(mask == 0)[0]  # Invert previous mask
    regions = [region for region, in ndimage.find_objects(mask)]
    msa3 = msa1.loc[:, regions]
    full2trim, i = {}, 0  # Full to trimmed MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            full2trim[j] = i
            i += 1

    # 2.3.2 Calculate gap scores of non-gap MSA
    scores = np.zeros(msa3.shape[1])
    for i, col in enumerate(msa3.iter_positions()):
        scores[i] = col.count('-')
    local_signal = ndimage.gaussian_filter1d(scores / len(msa3), sigma=l_sigma, mode='constant', cval=1)

    # 2.4 Make signals arrays (for nonlocal gap bias and trimming near gaps)
    gap_signals = np.zeros(msa1.shape)
    trim_signals = np.zeros(msa1.shape)
    gaps_array1 = np.array([[sym == '-' for sym in str(seq)] for seq in msa1])  # Sequences coded as boolean gap or non-gap

    # 2.5 Get trim slices
    trims1 = []
    for segment in segments:
        # 2.5.1 Unpack names
        region = segment['region']
        index = segment['index']
        slices = segment['slices']
        length = sum([s.stop-s.start for s in slices])

        # 2.5.2 Get local gap bias
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

        # 2.5.3 Get nonlocal gap bias
        gap_signal = gap_signals[index]
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        nonlocal_bias = gap_signal[start:stop].sum()

        # 2.5.4 Update trims, nonlocal bias, and trim signal
        if is_trimmed(length, segment['support'], local_bias, nonlocal_bias):
            trims1.extend(slices)
            propagate(start, stop, length, gap_signal, gaps_array1[index], nl_rate, nl_cutoff)  # Nonlocal bias
            propagate(start, stop, length, trim_signals[index], gaps_array1[index], tg_rate, tg_cutoff)  # Trim signal
    trims1 = sorted(trims1, key=lambda x: (x.start, x.start-x.stop))  # Start then -length

    # 2.6 Merge and invert slices
    # 2.6.1 Merge slices
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

    # 2.6.2 Invert slices
    trims3 = []
    start = 0
    for trim in trims2:
        trims3.append(slice(start, trim.start))
        start = trim.stop
    trims3.append(slice(start, msa3.shape[1]))

    # 2.7 Combine trims (segments and columns) to yield final alignment
    gaps_array2 = trim_signals > tg_threshold
    for i, gaps2 in enumerate(gaps_array2):
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            syms[region.start:region.stop] = (region.stop - region.start) * ['-']
    seqs = [skbio.Protein(''.join(syms), metadata=seq.metadata) for syms, seq in zip(syms_list, msa1)]
    msa4 = skbio.TabularMSA(seqs).loc[:, trims3]
    return msa4


def propagate(start, stop, length, signal, gaps, rate, cutoff):
    """Propagate an exponentially decreasing signal.

    The signal decays over gaps; however, the signal is only added to
    coordinates corresponding non-gap symbols.
    """
    max_k = -log(cutoff / length) / rate

    k = 0
    while start - k - 1 >= 0 and k < max_k:
        if not gaps[start - k - 1]:
            v = length * exp(-tc_rate * k)
            signal[start - k - 1] += v
        k += 1

    k = 0
    while stop + k <= len(signal) - 1 and k < max_k:
        if not gaps[stop + k]:
            v = length * exp(-tc_rate * k)
            signal[stop + k] += v
        k += 1


def get_segments(msa, region):
    """Return non-gap segments of a sequence in a region of an MSA."""
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


def is_trimmed(length, support, local_bias, nonlocal_bias):
    """Return boolean of whether to trim segment using logistic function."""
    b = 0
    w1, w2, w3, w4 = 1, -1, -1, 2
    threshold = 0.5

    x = b + w1*length ** 2 + w2*support + w3*local_bias + w4*nonlocal_bias
    p = 1 / (1 + exp(-x))
    return p > threshold
