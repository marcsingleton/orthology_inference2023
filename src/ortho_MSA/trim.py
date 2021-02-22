"""Trim columns and segments from an alignment."""

from math import exp, log

import numpy as np
import scipy.ndimage as ndimage
import skbio
from src.ortho_MSA.constants import constants


def trim_msa(msa):
    """Trim MSA by removing segments near indels in conserved regions and large insertions."""
    # 1 Calculate shared variables
    scores = np.zeros(msa.shape[1])
    for i, col in enumerate(msa.iter_positions()):
        scores[i] = col.count('-')
    gaps_array = np.array([[sym == '-' for sym in str(seq)] for seq in msa])

    # 2 Get trims (segments and columns)
    syms_list1 = trim_conserved(msa, scores, gaps_array)
    syms_list2, trims = trim_insertions(msa, scores, gaps_array)

    # 3.1 Extract slices from column trims
    slices = []
    for trim in trims:
        if trim['trimmed']:
            slices.extend(trim['slices'])
    slices = sorted(slices, key=lambda x: (x.start, x.start - x.stop))  # Start then -length

    # 3.2 Merge slices
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

    # 3.3 Invert slices
    slices_invert = []
    start = 0
    for s in slices_invert:
        slices_invert.append(slice(start, s.start))
        start = s.stop
    slices_invert.append(slice(start, msa.shape[1]))

    # 4 Combine trims (segments and columns) to yield final alignment
    seqs = []
    for seq, syms1, syms2 in zip(msa, syms_list1, syms_list2):
        syms = ['-' if sym1 != sym2 else sym1 for sym1, sym2 in zip(syms1, syms2)]  # Will only differ if one is converted to gap
        seqs.append(skbio.Protein(''.join(syms), metadata=seq.metadata))
    return skbio.TabularMSA(seqs).loc[:, slices_invert]


def trim_conserved(msa1, scores1, gaps_array1):
    """Trim MSA by removing segments near indels in conserved regions."""
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get conserved regions
    binary = ndimage.binary_closing(scores1 <= len(msa1) * constants['CON_FRAC'])
    mask = ndimage.label(binary, structure=constants['CON_CLOSE'] * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= constants['CON_MINLEN']]
    trim2full, i = {}, 0  # Trimmed to full MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            trim2full[i] = j
            i += 1

    msa2 = msa1.loc[:, regions]
    gaps_array2 = np.concatenate([gaps_array1[:, region] for region in regions], axis=1)

    # 2 Trim segments near gaps
    for i, seq in enumerate(msa2):
        # 2.1 Get gap segments
        gaps2 = gaps_array2[i]  # Sequence coded as boolean gap or non-gap
        trim_signal = np.zeros(len(seq))
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            length = region.stop - region.start
            propagate(region.start, region.stop, length, trim_signal, gaps2, constants['CON_RATE'])

        # 2.2 Trim non-gap segments
        gaps3 = trim_signal >= constants['CON_MINSIG']  # Sequence coded as boolean gap or non-gap after signal propagation
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps3)[0]):
            for j in range(region.start, region.stop):  # Iterate over positions to not delete segments between boundaries
                syms[trim2full[j]] = '-'
    return syms_list


def trim_insertions(msa1, scores1, gaps_array1, trimmed_dict=None):
    """Trim MSA by removing large insertions."""
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get gap regions
    mask = ndimage.label(len(msa1) - scores1 <= constants['GAP_NUM'])[0]
    regions = [region for region, in ndimage.find_objects(mask)]

    # 2 Get segments in regions
    segments = []
    for region in regions:
        segments.extend(get_segments(msa1, region))
    segments = sorted(segments, key=lambda x: sum([s.stop-s.start for s in x['slices']]), reverse=True)

    # 3 Make signal vector for local gap bias
    mask = ndimage.label(mask == 0)[0]  # Invert previous mask
    regions = [region for region, in ndimage.find_objects(mask)]
    full2trim, i = {}, 0  # Full to trimmed MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            full2trim[j] = i
            i += 1

    scores2 = np.concatenate([scores1[region] for region in regions])
    local_signal = ndimage.gaussian_filter1d(scores2 / len(msa1), sigma=constants['LOCAL_SIGMA'], mode='constant', cval=1)

    # 4 Make signals arrays for nonlocal gap bias and trims near gaps
    nonlocal_signals = np.zeros(msa1.shape)
    trim_signals = np.zeros(msa1.shape)

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
        nonlocal_signal = nonlocal_signals[index]
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        nonlocal_bias = nonlocal_signal[start:stop].sum()

        # 5.4 Update trims, nonlocal signal, and trim signal
        if trimmed_dict is None:
            trimmed = is_trimmed(length, segment['support'], local_bias, nonlocal_bias)
        else:
            trimmed = trimmed_dict.get((region.start, region.stop, index), None)
        trims.append({'region': region, 'index': index, 'slices': slices,
                      'trimmed': trimmed,
                      'length': length, 'support': support, 'local_bias': local_bias, 'nonlocal_bias': nonlocal_bias})
        if trimmed:
            propagate(start, stop, length, nonlocal_signal, gaps_array1[index], constants['NONLOCAL_RATE'])
            propagate(start, stop, length, trim_signals[index], gaps_array1[index], constants['GAP_RATE'])

    # 6 Trim segments
    gaps_array2 = trim_signals > constants['GAP_MINSIG']
    for i, gaps2 in enumerate(gaps_array2):
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps2)[0]):
            syms[region.start:region.stop] = (region.stop - region.start) * ['-']
    return syms_list, trims


def propagate(start, stop, length, signal, gaps, rate):
    """Propagate an exponentially decreasing signal.

    The signal decays over gaps; however, the signal is only added to
    coordinates corresponding non-gap symbols.
    """
    max_k = -log(1E-3 / length) / rate

    k = 0
    while start - k - 1 >= 0 and k < max_k:
        if not gaps[start - k - 1]:
            v = length * exp(-rate * k)
            signal[start - k - 1] += v
        k += 1

    k = 0
    while stop + k <= len(signal) - 1 and k < max_k:
        if not gaps[stop + k]:
            v = length * exp(-rate * k)
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
    weights = [constants['W0'], constants['W1'], constants['W2'], constants['W3'], constants['W4']]
    regressors = [1, length**2, support, local_bias, nonlocal_bias]
    x = sum([w*r for w, r in zip(weights, regressors)])
    p = 1 / (1 + exp(-x))
    return p > constants['THRESHOLD']
