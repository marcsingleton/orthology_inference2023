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

# LOGISTIC CLASSIFIER PARAMETERS
w0, w1, w2, w3, w4 = 0, 1, -1, -1, 2
threshold = 0.5


def trim_msa(msa1):
    """Trim MSA by removing large insertions and segments near indels."""
    scores = np.zeros(msa1.shape[1])
    for i, col in enumerate(msa1.iter_positions()):
        scores[i] = col.count('-')

    syms_list1 = trim_conserved(msa1, scores)
    syms_list2, trims = trim_insertions(msa1, scores)

    # Extract slices
    slices = []
    for trim in trims:
        if trim['trimmed']:
            slices.extend(trim['slices'])
    slices = sorted(slices, key=lambda x: (x.start, x.start - x.stop))  # Start then -length

    # Merge slices
    slices_merge = []
    if slices:
        start, stop = slices[0].start, slices[0].stop
        for s in slices:
            if s.start > stop:
                slices_merge.append(slice(start, stop))
                start = s.start
            if s.stop > stop:
                stop = s.stop
        slices_merge.append(slice(start, stop))

    # Invert slices
    slices_invert = []
    start = 0
    for s in slices_invert:
        slices_invert.append(slice(start, s.start))
        start = s.stop
    slices_invert.append(slice(start, msa1.shape[1]))

    # Combine trims (segments and columns) to yield final alignment
    seqs = []
    for seq, syms1, syms2 in zip(msa1, syms_list1, syms_list2):
        syms = ['-' if sym1 != sym2 else sym1 for sym1, sym2 in zip(syms1, syms2)]  # Will only differ if one is converted to gap
        seqs.append(skbio.Protein(''.join(syms), metadata=seq.metadata))
    msa2 = skbio.TabularMSA(seqs).loc[:, slices_invert]
    return msa2


def trim_conserved(msa1, scores):
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get conserved regions
    mask = ndimage.label(ndimage.binary_closing(scores / len(msa1) < c_fraction), structure=c_close * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= c_length]
    msa2 = msa1.loc[:, regions]
    trim2full, i = {}, 0  # Trimmed to full MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            trim2full[i] = j
            i += 1

    # 2 Trim segments near gaps
    for i, seq in enumerate(msa2):
        # 2.1 Get gap segments
        gaps1 = np.array([sym == '-' for sym in str(seq)])  # Sequence coded as boolean gap or non-gap
        trim_signal = np.zeros(len(gaps1))  # Deletion signals
        for region, in ndimage.find_objects(ndimage.label(gaps1)[0]):
            length = region.stop - region.start
            propagate(region.start, region.stop, length, trim_signal, gaps1, tc_rate, tc_cutoff)

        # 2.2 Trim non-gap segments
        gaps2 = trim_signal > tc_threshold  # Sequence coded as boolean gap or non-gap after signal propagation
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            for j in range(region.start, region.stop):  # Iterate over positions to not delete segments between boundaries
                syms[trim2full[j]] = '-'
    return syms_list


def trim_insertions(msa1, scores):
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get gap regions
    mask = ndimage.label(scores / len(msa1) > g_fraction)[0]
    regions = [region for region, in ndimage.find_objects(mask)]

    # 2 Get segments in regions
    segments = []
    for region in regions:
        segments.extend(get_segments(msa1, region))
    segments = sorted(segments, key=lambda x: sum([s.stop-s.start for s in x['slices']]), reverse=True)

    # 3 Get local gap bias
    # 3.1 Make MSA of non-gap regions
    mask = ndimage.label(mask == 0)[0]  # Invert previous mask
    regions = [region for region, in ndimage.find_objects(mask)]
    msa2 = msa1.loc[:, regions]
    full2trim, i = {}, 0  # Full to trimmed MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            full2trim[j] = i
            i += 1

    # 3.2 Calculate gap scores of non-gap MSA
    scores = np.zeros(msa2.shape[1])
    for i, col in enumerate(msa2.iter_positions()):
        scores[i] = col.count('-')
    local_signal = ndimage.gaussian_filter1d(scores / len(msa2), sigma=l_sigma, mode='constant', cval=1)

    # 4 Make signals arrays (for nonlocal gap bias and trimming near gaps)
    gap_signals = np.zeros(msa1.shape)
    trim_signals = np.zeros(msa1.shape)
    gaps_array1 = np.array([[sym == '-' for sym in str(seq)] for seq in msa1])  # Sequences coded as boolean gap or non-gap

    # 5 Get trim slices
    trims = []
    for segment in segments:
        # 5.1 Unpack names
        region = segment['region']
        index = segment['index']
        slices = segment['slices']
        support = segment['support']
        length = sum([s.stop-s.start for s in slices])

        # 5.2 Get local gap bias
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

        # 5.3 Get nonlocal gap bias
        gap_signal = gap_signals[index]
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        nonlocal_bias = gap_signal[start:stop].sum()

        # 5.4 Update trims, nonlocal bias, and trim signal
        trimmed = is_trimmed(length, segment['support'], local_bias, nonlocal_bias)
        trims.append({'region': region, 'index': index, 'slices': slices,
                      'trimmed': trimmed,
                      'length': length, 'support': support, 'local_bias': local_bias, 'nonlocal_bias': nonlocal_bias})
        if trimmed:
            propagate(start, stop, length, gap_signal, gaps_array1[index], nl_rate, nl_cutoff)  # Nonlocal bias
            propagate(start, stop, length, trim_signals[index], gaps_array1[index], tg_rate, tg_cutoff)  # Trim signal

    # 6 Trim segments
    gaps_array2 = trim_signals > tg_threshold
    for i, gaps2 in enumerate(gaps_array2):
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            syms[region.start:region.stop] = (region.stop - region.start) * ['-']
    return syms_list, trims


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
    x = w0 + w1*length**2 + w2*support + w3*local_bias + w4*nonlocal_bias
    p = 1 / (1 + exp(-x))
    return p > threshold
