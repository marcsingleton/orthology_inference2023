"""Functions for trimming columns from alignments."""

import scipy.ndimage as ndimage


def get_complement_slices(slices, start=0, stop=None):
    """Return slices of complement of input slices.

    Parameters
    ----------
    slices: list of slices
        Must be sorted and merged.
    start: int
        Start of interval from which complement slices are given.
    stop: int
        Stop of interval from which complement slices are given.

    Returns
    -------
    complement: list of slices
    """
    complement = []
    if slices:
        start0, stop0 = slices[0].start, slices[0].stop
        if start < start0:
            complement.append(slice(start, start0))
        for s in slices[1:]:
            complement.append(slice(stop0, s.start))
            stop0 = s.stop
        if stop is None or stop0 < stop:
            complement.append(slice(stop0, stop))
    else:
        complement.append(slice(start, stop))
    return complement


def get_merged_slices(slices):
    """Return slices where overlapping slices are merged.

    Parameters
    ----------
    slices: list of slices

    Returns
    -------
    merged: list of slices
    """
    merged = []
    if slices:
        slices = sorted(slices, key=lambda x: x.start)
        start0, stop0 = slices[0].start, slices[0].stop
        for s in slices[1:]:
            if s.start > stop0:
                merged.append(slice(start0, stop0))
                start0, stop0 = s.start, s.stop
            elif s.stop > stop0:
                stop0 = s.stop
        merged.append((slice(start0, stop0)))  # Append final slice
    return merged


def get_bound(profile, gradient, start, stop, sign):
    """Return column index with the largest weighted change in number of gaps in the range from start to stop.

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
    profile: 1-dimensional ndarray
        Weighted counts of gap symbols in MSA.
    gradient: 1-dimensional ndarray
        Gradient of posterior decoding. Weights the raw changes in gap number
        to choose the column with the sharpest inflection point.
    start: int
        Start of interval to test.
    stop: int
        Stop of interval to test, non-inclusive.
    sign: +1 or -1
        -1 for left bounds and +1 for right bounds.

    Returns
    -------
    j: int
    """
    deltas = []
    for i in range(start+1, stop):
        deltas.append((i, abs(gradient[i]*(profile[i] - profile[i-1]))))
    return max(deltas, key=lambda x: (x[1], sign*x[0]))[0]


def get_hull_slices(posterior, gradient, posterior_high, posterior_low, gradient_low):
    """Return slices corresponding to the hull that cover the initial region found in get_trim_slices.

    Parameters
    ----------
    posterior: 1-dimensional ndarray
        Posterior decoding.
    gradient: 1-dimensional ndarray
        Gradient of posterior decoding.
    posterior_high: float
        Cutoff for core region.
    posterior_low: float
        Cutoff for expanding margin outward based on posterior.
    gradient_low: float
        Cutoff for expanding margin outward based on gradient.
    """
    slices = []
    for s0, in ndimage.find_objects(ndimage.label(posterior >= posterior_high)[0]):
        start = s0.start  # start of left margin
        while start-1 >= 0 and (posterior[start-1] >= posterior_low or gradient[start-1] >= gradient_low):
            start -= 1
        stop = s0.stop - 1  # stop of right margin
        while stop+1 < len(posterior) and (posterior[stop+1] >= posterior_low or gradient[stop+1] <= -gradient_low):
            stop += 1
        slices.append(slice(start, stop+1))
    get_merged_slices(slices)

    return slices


def get_trim_slices(profile, posterior, gradient, posterior_high, posterior_low, gradient_high, gradient_low):
    """Return slices corresponding to trimmed regions.

    This algorithm is designed to cleanly trim alignments at points where the
    number of gaps changes rapidly. It does this by first identifying core
    regions which exceed the posterior_high cutoff. It identifies margins on
    the left and right sides of this core region. The extent of the margins is
    determined by the posterior_low, gradient_high, and gradient_low
    cutoffs. When expanding the margin outward, columns are added if the
    posterior or gradient is higher than posterior_low and gradient_low,
    respectively. When expanding inward, columns are added if the gradient is
    higher than gradient_high. These rules ensure that the margins of the core
    region fully capture the shoulders of posterior peaks.

    Specific columns from these left and right margins are then chosen using
    the get_bound function.

    Parameters
    ----------
    profile: 1-dimensional ndarray
        Weighted counts of gap symbols in MSA.
    posterior: 1-dimensional ndarray
        Posterior decoding.
    gradient: 1-dimensional ndarray
        Gradient of posterior decoding.
    posterior_high: float
        Cutoff for core region.
    posterior_low: float
        Cutoff for expanding margin outward based on posterior.
    gradient_high: float
        Cutoff for expanding margin inward based on gradient.
    gradient_low: float
        Cutoff for expanding margin outward based on gradient.
    """
    slices = []
    for s0, in ndimage.find_objects(ndimage.label(posterior >= posterior_high)[0]):
        lstart = s0.start  # start of left margin
        while lstart-1 >= 0 and (posterior[lstart-1] >= posterior_low or gradient[lstart-1] >= gradient_low):
            lstart -= 1
        lstop = s0.start  # stop of left margin
        while lstop+1 < len(posterior) and gradient[lstop+1] >= gradient_high:
            lstop += 1
        if lstart < lstop:
            start = get_bound(profile, gradient, lstart, lstop+1, -1)
        else:
            start = s0.start

        rstart = s0.stop - 1  # start of right margin
        while rstart-1 >= 0 and gradient[rstart-1] <= -gradient_high:
            rstart -= 1
        rstop = s0.stop - 1  # stop of right margin
        while rstop+1 < len(posterior) and (posterior[rstop+1] >= posterior_low or gradient[rstop+1] <= -gradient_low):
            rstop += 1
        if rstop > s0.stop:
            stop = get_bound(profile, gradient, rstart, rstop+1, 1)
        else:
            stop = s0.stop

        if start < stop:
            s = slice(start, stop)
        else:
            s = s0
        slices.append(s)
    slices = get_merged_slices(slices)

    return slices
