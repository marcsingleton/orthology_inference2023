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


def get_bound(msa, gradient, start, stop, sign):
    """Return column with largest change in number of gaps from previous.

    Gaps are calculated from the range of lower to upper, inclusive. Ties are
    broken towards the most "exterior" column. (The sign argument, which is
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
        deltas.append((j, abs(gradient[j])*(count - count0)))
        count0 = count
    return max(deltas, key=lambda x: (sign*x[1], -sign*x[0]))[0]


def get_slices(msa, posterior, gradient):
    """Return slices corresponding to trimmed regions.

    The first and last column of the HIGH cutoff are included so at least one
    delta is calculated. Since the deltas are calculated as current - previous,
    the stop bound corresponds to the column with (presumably) a small number
    of gaps and will be correctly excluded from the slice.
    """
    slices = []
    for region, in ndimage.find_objects(ndimage.label(posterior >= POST_HIGH)[0]):
        lstart = region.start  # start of left margin
        while lstart-1 >= 0 and (posterior[lstart-1] >= POST_LOW or gradient[lstart-1] >= GRAD_LOW):
            lstart -= 1
        lstop = region.start  # stop of left margin
        while lstop+1 < len(posterior) and gradient[lstop+1] >= GRAD_HIGH:
            lstop += 1
        if lstart < lstop:
            start = get_bound(msa, gradient, lstart, lstop, 1)
        else:
            start = region.start

        rstart = region.stop - 1  # start of right margin
        while rstart-1 >= 0 and gradient[rstart-1] <= -GRAD_HIGH:
            rstart -= 1
        rstop = region.stop - 1  # stop of right margin
        while rstop+1 < len(posterior) and (posterior[rstop+1] >= POST_LOW or gradient[rstop+1] <= -GRAD_LOW):
            rstop += 1
        if rstop > region.stop:
            stop = get_bound(msa, gradient, rstart, rstop, -1)
        else:
            stop = region.stop

        if start < stop:
            s = slice(start, stop)
        else:
            s = region
        slices.append(s)

    return slices


POST_HIGH = 0.75
POST_LOW = 0.5
GRAD_HIGH = 0.02
GRAD_LOW = 0.001
