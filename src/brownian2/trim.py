"""Trim columns from alignment using decoded states from HMM."""

import scipy.ndimage as ndimage


def trim_terminals(msa):
    """Return MSA with unaligned terminal regions removed."""
    idx = 0
    for j in range(len(msa[0][1])):
        for i in range(len(msa)):
            sym = msa[i][1][j]
            if sym == '.' or sym.islower():
                break
        else:
            idx = j
            break  # if no break exit
    msa = [(header, seq[idx:]) for header, seq in msa]

    idx = len(msa[0][1])
    for j in range(len(msa[0][1]), 0, -1):
        for i in range(len(msa)):
            sym = msa[i][1][j - 1]
            if sym == '.' or sym.islower():
                break
        else:
            idx = j
            break  # if no break exit
    msa = [(header, seq[:idx]) for header, seq in msa]

    return msa


def get_bound(msa, start, stop, sign):
    """Return column with largest change in number of gaps from previous.

    Gaps are calculated from the range of lower to upper, inclusive. Ties are
    broken towards the most "interior" column. (The sign argument, which is
    either +1 or -1, correctly controls the ordering of both the delta and
    index for both the start and stop bounds.)
    """
    gaps = []
    for j in range(start, stop+1):
        count = 0
        for _, seq in msa:
            if seq[j] in ['-', '.']:
                count += 1
        gaps.append((j, count))

    deltas, count0 = [], gaps[0][1]
    for j, count in gaps[1:]:
        deltas.append((j, count - count0))
        count0 = count
    return max(deltas, key=lambda x: (sign*x[1], sign*x[0]))[0]


def get_slices(msa, posterior, gradient):
    """Return slices corresponding to trimmed regions.

    The first and last column of the HIGH cutoff are included so at least one
    delta is calculated. Since the deltas are calculated as current - previous,
    the stop bound corresponds to the column with (presumably) a small number
    of gaps and will be correctly excluded from the slice.
    """
    slices = []
    for region, in ndimage.find_objects(ndimage.label(posterior >= HIGH)[0]):
        start = region.start  # start of left margin
        while start-1 >= 0 and posterior[start-1] >= LOW and gradient[start-1] >= GRADIENT:
            start -= 1
        if start < region.start:
            start = get_bound(msa, start, region.start, 1)

        stop = region.stop  # stop of right margin
        while stop+1 < len(posterior) and posterior[stop+1] >= LOW and gradient[stop+1] <= -GRADIENT:
            stop += 1
        if stop > region.stop:
            stop = get_bound(msa, region.stop - 1, stop, -1)

        slices.append(slice(start, stop))

    return slices


HIGH = 0.6
LOW = 0.01
GRADIENT = 0.005
