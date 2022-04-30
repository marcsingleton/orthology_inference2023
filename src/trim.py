"""Functions for trimming columns from alignments."""

import scipy.ndimage as ndimage


def trim_terminals(msa):
    """Return MSA with unaligned terminal regions removed.

    Unaligned regions are created by HMMer a symbol in the sequence does not
    match any of the consensus columns columns in the profile HMM. These
    regions are indicated with lowercase symbols for . for gaps.

    Parameters
    ----------
    msa: list of (header, seq)

    Returns
    -------
    msa: list of (header, seq)
    """
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
    """Return column index with largest weighted change in number of gaps in range from start to stop.

    Gaps are calculated from the range of lower to upper, inclusive. Ties are
    broken towards the most "exterior" column via the sign argument, which is
    either +1 or -1, and correctly controls the ordering of both the delta and
    index for both the start and stop bounds.)

    Parameters
    ----------
    msa: list of (header, seq)
    gradient: 1-dimenionsal ndarray
        Gradient of posterior decoding. Weights the raw changes in gap number
        to choose the column with the sharpest inflection point.
    start: int
        Start of interval to test.
    stop: int
        Stop of interval to test.
    sign: +1 or -1
        +1 for left bounds and -1 for right bounds.

    Returns
    -------
    j: int
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


def get_slices(msa, posterior, gradient, posterior_high, posterior_low, gradient_high, gradient_low):
    """Return slices corresponding to trimmed regions.

    The first and last column of the gradient_high threshold are included so at
    least one delta is calculated. Since the deltas are calculated as
    current - previous, the stop bound corresponds to the column with
    (presumably) a small number of gaps and will be correctly excluded from the
    slice.

    Parameters
    ----------
    msa: list of (header, seq)
    posterior: 1-dimenionsal ndarray
        Posterior decoding.
    gradient: 1-dimenionsal ndarray
        Gradient of posterior decoding.
    posterior_high: float
        Threshold for core region.
    posterior_low: float
        Threshold for expanding margin outward based on posterior.
    gradient_high: float
        Threshold for expanding margin inward based on gradient.
    gradient_low: float
        Threshold for expanding margin outward based on gradient.
    """
    slices = []
    for region, in ndimage.find_objects(ndimage.label(posterior >= posterior_high)[0]):
        lstart = region.start  # start of left margin
        while lstart-1 >= 0 and (posterior[lstart-1] >= posterior_low or gradient[lstart - 1] >= gradient_low):
            lstart -= 1
        lstop = region.start  # stop of left margin
        while lstop+1 < len(posterior) and gradient[lstop+1] >= gradient_high:
            lstop += 1
        if lstart < lstop:
            start = get_bound(msa, gradient, lstart, lstop, 1)
        else:
            start = region.start

        rstart = region.stop - 1  # start of right margin
        while rstart-1 >= 0 and gradient[rstart-1] <= -gradient_high:
            rstart -= 1
        rstop = region.stop - 1  # stop of right margin
        while rstop+1 < len(posterior) and (posterior[rstop+1] >= posterior_low or gradient[rstop + 1] <= -gradient_low):
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
