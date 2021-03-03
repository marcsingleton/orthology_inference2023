"""Trim columns and segments from an alignment."""

from itertools import combinations
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
    slices_invert, i = [], 0
    for s in slices_invert:
        slices_invert.append(slice(i, s.start))
        i = s.stop
    slices_invert.append(slice(i, msa.shape[1]))

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
    mask = ndimage.label(binary, structure=constants['CON_WINDOW'] * [1])[0]
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
    gap_propensity = ndimage.gaussian_filter1d(scores2 / len(msa1), sigma=constants['GP_SIGMA'], mode='constant', cval=0)

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
        propagate(start, stop, length, indel_signals1[index], gaps_array1[index], constants['INDEL1_RATE'])

    indel_signals2 = np.zeros(msa1.shape)
    binary = ndimage.binary_closing(scores1 <= len(msa1) * constants['CON_FRAC'])
    mask = ndimage.label(binary, structure=constants['CON_WINDOW'] * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= constants['CON_MINLEN']]
    for region1 in regions:
        gaps_array2 = gaps_array1[:, region1]
        for i, gaps2 in enumerate(gaps_array2):
            for region2, in ndimage.find_objects(ndimage.label(gaps2)[0]):
                propagate(region1.start + region2.start, region1.start + region2.stop,
                          region2.stop - region2.start, indel_signals2[i], gaps_array1[i], constants['INDEL2_RATE'])

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
        w = (constants['GD_WINDOW'] - 1) // 2
        n = 2**len(msa1) * (1 - ((2**len(msa1) - 1)/(2**len(msa1)))**constants['GD_WINDOW'])
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
        if trimmed_dict is None:
            trimmed = is_trimmed(length, segment['support'], gp, gd, ib1, ib2)
        else:
            trimmed = trimmed_dict.get((region.start, region.stop, index), None)
        trims.append({'region': region, 'index': index, 'slices': slices, 'trimmed': trimmed,
                      'length': length, 'support': support,
                      'gap_propensity': gp, 'gap_diversity': gd,
                      'indel_bias1': ib1, 'indel_bias2': ib2})
        if trimmed:
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
    truncate_k = int(-log(1E-3 / length) / rate) + 1

    max_k = min(truncate_k, start)  # Number of positions to propagate over
    s = [length * exp(rate * (k-max_k+1)) for k in range(max_k)]
    signal[start-max_k:start] += s * ~gaps[start-max_k:start]

    max_k = min(truncate_k, len(signal) - stop)  # Number of positions to propagate over
    s = [length * exp(-rate * k) for k in range(max_k)]
    signal[stop:stop+max_k] += s * ~gaps[stop:stop+max_k]


def get_segments(msa, region):
    """Return non-gap segments of a sequence in a region of an MSA."""
    block = msa.loc[:, region]
    segments = {i: {'region': region, 'index': i, 'slices': [], 'support': 0} for i in range(len(block))}
    starts = {i: None for i in range(len(block))}  # Slice starts
    for j, col in enumerate(block.iter_positions()):
        support = 0
        syms = [sym for sym in str(col) if sym != '-']
        for sym1, sym2 in combinations([sym for sym in str(col) if sym != '-'], 2):
            support += constants['MATRIX'][(sym1, sym2)]
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


def is_trimmed(length, support, gap_propensity, gap_diversity, indel_bias1, indel_bias2):
    """Return boolean of whether to trim segment using logistic function."""
    weights = [constants['W0'], constants['W1'], constants['W2'], constants['W3'], constants['W4'], constants['W5'], constants['W6']]
    regressors = [1, length, support, gap_propensity, gap_diversity, indel_bias1, indel_bias2]
    x = sum([w*r for w, r in zip(weights, regressors)])
    p = 1 / (1 + exp(-x))
    return p > constants['THRESHOLD']
