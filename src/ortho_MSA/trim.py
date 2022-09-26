"""Functions for trimming columns from alignments."""

import scipy.ndimage as ndimage


def get_bound(msa, gradient, start, stop, sign):
    """Return column index with largest weighted change in number of gaps in the range from start to stop.

    To select the column where the number of gaps changes most rapidly, the
    gradient data is multiplied by the raw gap deltas. This weighting should
    select the column with a large gap delta that is also supported by the
    posterior decoding, i.e. its confidence is also rapidly changing at that
    column.

    Ties are broken towards the most "exterior" column via the sign argument,
    which is either +1 or -1 and correctly controls the ordering of both the
    delta and index for both the start and stop bounds.

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

    This algorithm is designed to cleanly trim alignments at points where the
    number of gaps changes rapidly. It does this by first identifying core
    regions which exceed the posterior_high threshold. It identifies margins on
    the left and right sides of this core region. The extent of the margins is
    determined by the posterior_low, gradient_high, and gradient_low
    thresholds. When expanding the margin outward, columns are added if the
    posterior or gradient is higher than posterior_low and gradient_low,
    respectively. When expanding inward, columns are added if the gradient is
    higher than gradient_high. These rules ensure that the margins of the core
    region fully capture the shoulders of posterior peaks.

    Specific columns from these left and right margins are then chosen using
    the get_bound function.

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
    slices = sorted(slices, key=lambda x: x.start)

    # Merge slices
    if slices:
        merged = []
        start0, stop0 = slices[0].start, slices[0].stop
        for s in slices[1:]:
            if s.start > stop0:
                merged.append(slice(start0, stop0))
                start0, stop0 = s.start, s.stop
            elif s.stop > stop0:
                stop0 = s.stop
        merged.append((slice(start0, stop0)))  # Append final slice

    return slices
