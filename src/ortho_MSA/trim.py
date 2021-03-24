"""Trim columns and segments from an alignment."""

from itertools import combinations
from math import exp, log

import numpy as np
import scipy.ndimage as ndimage


def trim_conserved(msa1, scores1, gaps_array1,
                   con_frac, con_window, con_minlen, con_rate, con_minsig):
    """Trim MSA by removing segments near indels in conserved regions.

    Parameters
    ----------
    msa1: TabularMSA (skbio)
    scores1: list or 1D array
        1D iterable where each element counts the number of gaps in each column
        of msa1.
    gaps_array1: 2D array
        2D array with same shape as msa1 where each position is coded as 1
        (gap) or 0 (non-gap).
    con_frac: float
        Float between 0 and 1 setting the maximum gap fraction in conserved
        regions.
    con_window: int
        Size of window for closing operation which combines adjacent conserved
        regions. This is intended prevent the subsequent length filtering from
        removing conserved regions which are interrupted by small gaps.
    con_minlen: int
        Minimum size of conserved region.
    con_rate: float
        Decay rate of trim signal.
    con_minsig: float
        Minimum signal to "trim" a position into a gap.

    Returns
    -------
    syms_list: list of lists
        List of symbols (as list) for each position in msa1. The first and
        second indices correspond to rows and columns, respectively. These
        symbols differ from msa1 in that short non-gap segments near large gaps
        are replaced with gap symbols.
    """
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get conserved regions
    binary = ndimage.binary_closing(scores1 <= len(msa1) * con_frac)
    mask = ndimage.label(binary, structure=con_window * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= con_minlen]
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
            propagate(region.start, region.stop, length, trim_signal, gaps2, con_rate)

        # 2.2 Trim non-gap segments
        gaps3 = trim_signal >= con_minsig  # Sequence coded as boolean gap or non-gap after signal propagation
        syms = syms_list[i]
        for region, in ndimage.find_objects(ndimage.label(gaps3)[0]):
            for j in range(region.start, region.stop):  # Iterate over positions to not delete segments between boundaries
                syms[trim2full[j]] = '-'
    return syms_list


def trim_insertions(msa1, scores1, gaps_array1,
                    gap_num, gap_rate, gap_minsig,
                    nongap_frac, nongap_window, nongap_minlen,
                    gp_sigma, gd_window, indel1_rate, indel2_rate,
                    weights, threshold,
                    matrix):
    """Trim MSA by removing large insertions.

    nongap_frac, nongap_window, and nongap_minlen are defined as in trim_conserved.
    weights, threshold are defined in is_trimmed.
    matrix is defined in get_segments.

    Parameters
    ----------
    msa1: TabularMSA (skbio)
    scores1: list or 1D array
        1D iterable where each element counts the number of gaps in each column
        of msa1.
    gaps_array1: 2D array
        2D array with same shape as msa1 where each position is coded as 1
        (gap) or 0 (non-gap).
    gap_num: int
        The maximum number of non-gap symbols in a gap region.
    gp_sigma: float
        The smoothing radius for calculating gap propensity. Gap propensity is
        defined as the local fraction of gaps in the alignment.
    gd_window: int
        Size of window for calculating gap diversity. Gap diversity is defined
        as the local number of unique gap patterns divided by the expected
        number of unique gap patterns.
    indel1_rate: float
        Decay rate of insertions in gap regions.
    indel2_rate: float
        Decay rate of deletions in conserved regions.
    gap_rate: float
        Decay rate of trim signal.
    gap_minsig: float
        Minimum signal to "trim" a position into a gap.
    """
    syms_list = [list(str(seq)) for seq in msa1]   # List of symbols representing alignment with trimmed segments

    # 1 Get gap regions
    mask = ndimage.label(len(msa1) - scores1 <= gap_num)[0]
    regions = [region for region, in ndimage.find_objects(mask)]

    # 2 Get segments in regions
    segments = []
    for region in regions:
        segments.extend(get_segments(msa1, region, matrix))

    # 3.1 Make arrays for gap metrics
    mask = ndimage.label(mask == 0)[0]  # Invert previous mask
    regions = [region for region, in ndimage.find_objects(mask)]
    full2trim, i = {}, 0  # Full to trimmed MSA coordinates
    for region in regions:
        for j in range(region.start, region.stop):
            full2trim[j] = i
            i += 1

    # 3.2 Gap propensity
    scores2 = np.concatenate([scores1[region] for region in regions])
    gap_propensity = ndimage.gaussian_filter1d(scores2 / len(msa1), sigma=gp_sigma, mode='constant', cval=0)

    # 3.3 Gap diversity
    msa2 = msa1.loc[:, regions]
    gap_diversity = []
    for col in msa2.iter_positions():
        gap_diversity.append(frozenset([i for i, sym in enumerate(str(col)) if sym == '-']))

    # 4 Make arrays for indel and trim signals
    trim_signals = np.zeros(msa1.shape)
    indel_signals1 = np.zeros(msa1.shape)
    for segment in segments:
        index = segment['index']
        slices = segment['slices']
        length = sum([s.stop-s.start for s in slices])
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        propagate(start, stop, length, indel_signals1[index], gaps_array1[index], indel1_rate)

    indel_signals2 = np.zeros(msa1.shape)
    binary = ndimage.binary_closing(scores1 <= len(msa1) * nongap_frac)
    mask = ndimage.label(binary, structure=nongap_window * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= nongap_minlen]
    for region1 in regions:
        gaps_array2 = gaps_array1[:, region1]
        for i, gaps2 in enumerate(gaps_array2):
            for region2, in ndimage.find_objects(ndimage.label(gaps2)[0]):
                propagate(region1.start + region2.start, region1.start + region2.stop,
                          region2.stop - region2.start, indel_signals2[i], gaps_array1[i], indel2_rate)

    # 5 Get trim slices
    trims = []
    for segment in segments:
        # 5.1 Unpack names
        region = segment['region']
        index = segment['index']
        slices = segment['slices']
        support = segment['support']
        length = sum([s.stop-s.start for s in slices])

        # 5.2 Get gap propensity and diversity
        w = int((gd_window - 1) // 2)
        n = 2**len(msa1) * (1 - ((2**len(msa1) - 1)/(2**len(msa1)))**gd_window)
        if region.start == 0:
            gp1 = 0
            gd1 = set()
        else:
            idx1 = full2trim[region.start-1]
            gp1 = gap_propensity[idx1]
            gd1 = set(gap_diversity[idx1-w:idx1+1])
        if region.stop == msa1.shape[1]:
            gp2 = 0
            gd2 = set()
        else:
            idx2 = full2trim[region.stop]
            gp2 = gap_propensity[idx2]
            gd2 = set(gap_diversity[idx2-w:idx2+1])
        gp = (gp1 + gp2) / 2
        gd = len(gd1 | gd2) / n

        # 5.3 Get indel biases
        start, stop = min([s.start for s in slices]), max([s.stop for s in slices])
        ib1 = indel_signals1[index, start:stop].sum()
        ib2 = indel_signals2[index, start:stop].sum()

        # 5.4 Update trims, indel signal, and trim signal
        trimmed = is_trimmed(length, segment['support'], gp, gd, ib1, ib2,
                             weights, threshold)
        trims.append({'region': region, 'index': index, 'slices': slices, 'trimmed': trimmed,
                      'length': length, 'support': support,
                      'gap_propensity': gp, 'gap_diversity': gd,
                      'indel_bias1': ib1, 'indel_bias2': ib2})
        if trimmed:
            propagate(start, stop, length, trim_signals[index], gaps_array1[index], gap_rate)
            syms = syms_list[index]
            for s in slices:
                syms[s.start:s.stop] = (s.stop - s.start) * ['-']

    # 6 Trim segments
    gaps_array2 = trim_signals > gap_minsig
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
    truncate_k = int(-log(1E-3 / length) / rate) + 1

    max_k = min(truncate_k, start)  # Number of positions to propagate over
    s = [length * exp(rate * (k-max_k+1)) for k in range(max_k)]
    signal[start-max_k:start] += s * ~gaps[start-max_k:start]

    max_k = min(truncate_k, len(signal) - stop)  # Number of positions to propagate over
    s = [length * exp(-rate * k) for k in range(max_k)]
    signal[stop:stop+max_k] += s * ~gaps[stop:stop+max_k]


def get_segments(msa, region, matrix):
    """Return non-gap segments of a sequence in a region of an MSA.

    Parameters
    ----------
    msa: TabularMSA (skbio)
    region: slice
        Slice selecting the region within msa.
    matrix: dict
        Dictionary of substitution matrix used to calculate support. Symbol
        pairs are keyed as tuples in both possible orders.

    Returns
    -------
    segments: list of dicts
        Each segment is a dictionary storing information about the non-gap
        symbols of each sequence in the region of msa. Specifically, it stores
        the following parameters:
            region: slice selecting the region within msa.
            index: index selecting the sequence within msa.
            support: sum of sum-of-pairs scores across each non-gap symbol in
                the sequence.
            slices: list of slices selecting each contiguous stretch of non-gap
                symbols that composes the segment.
    """
    block = msa.loc[:, region]
    segments = {i: {'region': region, 'index': i, 'slices': [], 'support': 0} for i in range(len(block))}
    starts = {i: None for i in range(len(block))}  # Slice starts
    for j, col in enumerate(block.iter_positions()):
        support = 0
        syms = [sym for sym in str(col) if sym != '-']
        for sym1, sym2 in combinations([sym for sym in str(col) if sym != '-'], 2):
            support += matrix[(sym1, sym2)]
        support *= 2/(max(2, len(syms))-1)  # Scale score to number of non-gap symbols
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


def is_trimmed(length, support, gap_propensity, gap_diversity, indel_bias1, indel_bias2,
               weights, threshold):
    """Return boolean of whether to trim segment using logistic function."""
    keys = ['bias', 'length', 'support', 'gap_propensity', 'gap_diversity', 'indel_bias1', 'indel_bias2']
    regressors = [1, length, support, gap_propensity, gap_diversity, indel_bias1, indel_bias2]
    x = sum([weights[key]*r for key, r in zip(keys, regressors)])
    p = 1 / (1 + exp(-x))
    return p > threshold
