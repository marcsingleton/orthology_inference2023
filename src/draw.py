"""Drawing functions for specialized data."""

from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


def _get_aspect(block_cols, ROWS, COLS, get_dims):
    length, height = get_dims(block_cols, ROWS, COLS)
    return length / height


def _get_block_cols(ROWS, COLS, ASPECT, get_dims):
    # Use binary search to find interval containing optimal block_cols
    # If an interval endpoint exactly equals the ratio, choose that endpoint
    # Otherwise choose the interval where the deviation changes sign between the two endpoints
    interval = (1, COLS)
    while interval[1] - interval[0] > 1:
        i0, i2 = interval
        i1 = (i0 + i2) // 2
        deltas = [_get_aspect(i, ROWS, COLS, get_dims) - ASPECT for i in (i0, i1, i2)]

        if any([delta == 0 for delta in deltas]):
            _, i = min(zip(deltas, (i0, i1, i2)), key=lambda x: abs(x[0]))
            interval = (i, i)  # Make degenerate interval to fit with min
            break
        elif deltas[0] * deltas[1] < 0:
            interval = (i0, i1)
        elif deltas[1] * deltas[2] < 0:
            interval = (i1, i2)
        else:
            interval = (i0, i2)
            break
    block_cols = min(range(interval[0], interval[1] + 1), key=lambda x: abs(_get_aspect(x, ROWS, COLS, get_dims) - ASPECT))  # Choose value that minimizes difference

    # Ensure last block is at least 50% of block_cols
    if COLS % block_cols < 0.5 * block_cols:  # Guarantees at least two blocks
        blocks_im = COLS // block_cols  # Total blocks minus 1
        block_cols += ceil((COLS % block_cols) / blocks_im)  # Distribute excess to other blocks
    return block_cols


def draw_msa(msa,
             aspect=2.5, hspace=25, sym_length=7, sym_height=7,
             block_cols=None, sym2color=None, gap2color=None):
    """Draw alignment as PNG.

    Parameters
    ----------
    msa: list of strings
    aspect: float
        Aspect ratio (length:height) of image. The function will calculate the
        number of columns in each block that yields an image that best matches
        this aspect ratio.
    hspace: int
        Number of pixels between blocks.
    sym_length: int
        Number of pixels in length of symbol rectangle.
    sym_height: int
        Number of pixels in height of symbol rectangle.
    block_cols: int
        Number of columns in each block. Will override ratio if is not None.
    sym2color: dict
        Mapping of symbols to color hex codes.
    gap2color: dict
        Mapping of symbols to color hex codes of center dot. Intended for
        distinguishing aligned and unaligned gaps produced by HMMer.

    Returns
    -------
        im: ndarray
            3 channel array with dimensions (height, length, channels)
    """
    # Define functions and globals
    ROWS, COLS = len(msa), len(msa[0])
    ASPECT = aspect

    def get_dims(block_cols, ROWS, COLS):
        im_length = sym_length * block_cols  # Length of final image
        block_num = COLS // block_cols - (1 if COLS % block_cols == 0 else 0)  # Number of blocks in addition to the first
        im_height = (sym_height * ROWS + hspace) * block_num + sym_height * ROWS  # Height of final image
        return im_length, im_height

    # Set options
    if block_cols is None:
        block_cols = _get_block_cols(ROWS, COLS, ASPECT, get_dims)
    if sym2color is None:
        sym2color = default_sym2color
    if gap2color is None:
        gap2color = default_gap2color

    # Instantiate array and fill with values
    im_length, im_height = get_dims(block_cols, ROWS, COLS)
    im = np.full((im_height, im_length, 3), 255, dtype='uint8')
    for i, seq in enumerate(msa):
        for j, sym in enumerate(seq):
            # Position of symbol rectangle
            block = j // block_cols
            y = (sym_height * ROWS + hspace) * block + sym_height * i
            x = j % block_cols * sym_length

            # Create RGB tuple
            try:
                color = sym2color[sym]
            except KeyError as error:
                print(f"Warning: Symbol {error} not in dictionary. Using color for symbol 'X' in its place.")
                color = sym2color['X']
            rgb = [int(color[i:i+2], 16) for i in (0, 2, 4)]

            # Fill slice with color
            im[y:y+sym_height, x:x+sym_length] = rgb
            if sym in gap2color:
                color = gap2color[sym]
                rgb = [int(color[i:i+2], 16) for i in (0, 2, 4)]
                y1, y2 = y + ceil((sym_height - 2) / 2), y + ceil((sym_height + 1) / 2)
                x1, x2 = x + ceil((sym_length - 2) / 2), x + ceil((sym_length + 1) / 2)
                im[y1:y2, x1:x2] = rgb
    return im


def plot_msa_data(msa, data, figsize=(15, 6),
                  msa_labels=None, msa_labelsize=6, msa_ticklength=0, msa_tickwidth=0.5, msa_tickpad=1,
                  x_start=0, x_labelsize=6, y_labelsize=6,
                  height_ratio=1, hspace=0.5, sym_length=7, sym_height=7,
                  data_min=None, data_max=None,
                  msa_legend=False, data_labels=None, legend_kwargs=None,
                  block_cols=None, sym2color=None, gap2color=None):
    """Plot MSA with associated positional data as matplotlib figure.

    Parameters
    ----------
    msa: list of strings
    data: list of lists or two-dimensional ndarray
        Data series to be plotted where outer list is the series and the inner
        list is a value associated with each column of the msa.
    figsize: 2-tuple of floats
        Figure dimension (width, height) in inches.
    msa_labels: list of strings
        Strings of labels of sequences in MSA in the same order as in msa.
        len(msa) must match len(msa_labels). Empty strings are not displayed.
    msa_labelsize: float
        Font size of MSA labels.
    msa_ticklength: float
        Length of MSA ticks.
    msa_tickwidth: float
        Width of MSA ticks.
    msa_tickpad: float
        Padding between MSA ticks and labels.
    x_start: int
        Starting value of x-axis.
    x_labelsize: float
        Font size of x-axis labels.
    y_labelsize: float
        Font size of y-axis labels.
    height_ratio: float
        Height of data axes as fraction of height of block.
    hspace: float
        Padding between each data axes and subsequent block of MSA as fraction
        of height of block.
    sym_length: int
        Number of pixels in length of the rectangles for each symbol.
    sym_height: int
        Number of pixels in height of the rectangles for each symbol.
    data_min: float
        Minimum of y-axis across all data axes.
    data_max: float
        Maximum of y-axis across all data axes.
    msa_legend: bool
        True if MSA legend is drawn.
    data_labels: list of strings
        Labels for data series. If not None, will draw legend. Length does not
        need to match the number of series in data.
    legend_kwargs: dict
        Additional kwargs passed to fig.legend call.
    block_cols: int
        The number of columns in each block. Overrides determination of optimal
        number if not None.
    sym2color: dict
        Mapping of symbols to colors given as hex strings. These values
        determine the background color of the rectangles.
    gap2color: dict
        Mapping of symbols to colors given as hex strings. These values
        determine the foreground color of the dash in gap symbols.

    Returns
    -------
    fig: Figure (matplotlib)
    """
    # Define functions and globals
    ROWS, COLS = len(msa), len(msa[0])
    ASPECT = figsize[0] / figsize[1]

    def get_dims(block_cols, ROWS, COLS):
        plot_length = block_cols  # Length of final plot
        block_num = COLS // block_cols - (1 if COLS % block_cols == 0 else 0)  # Number of blocks in addition to the first
        plot_height = (1 + height_ratio + 2*hspace) * ROWS * block_num + (1 + height_ratio + hspace) * ROWS  # Height of final image
        return plot_length, plot_height

    # Set options
    if msa_labels is None:
        msa_labels = ROWS * ['']
    if len(msa_labels) != ROWS:
        raise RuntimeError('len(msa_labels) does not match len(msa)')
    if legend_kwargs is None:
        legend_kwargs = {}
    if block_cols is None:
        block_cols = _get_block_cols(ROWS, COLS, ASPECT, get_dims)
    if sym2color is None:
        sym2color = default_sym2color
    if gap2color is None:
        gap2color = default_gap2color
    if isinstance(data, list):
        data = np.array(data)
    if data.ndim == 1:
        data = np.expand_dims(data, axis=0)
    if data_min is None:
        data_min = data.min() - 0.05 * (data.max() - data.min())
    if data_max is None:
        data_max = data.max() + 0.05 * (data.max() - data.min())
    if data_min == data_max:
        data_min -= 0.5
        data_max += 0.5
    block_num = COLS // block_cols + (1 if COLS % block_cols > 0 else 0)  # Number of blocks
    block_rows = len(msa)

    # Draw axes
    height_ratios = []
    for i in range(3*block_num):
        r = i % 3
        if r == 0:
            h = 1
        elif r == 1:
            h = height_ratio
        else:
            h = hspace
        height_ratios.append(h)

    im = draw_msa(msa, block_cols=len(msa[0]), sym_length=sym_length, sym_height=sym_height, sym2color=sym2color, gap2color=gap2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3*block_num, 1, figure=fig, height_ratios=height_ratios)
    for i in range(block_num):
        msa_ax = fig.add_subplot(gs[3*i, :])
        data_ax = fig.add_subplot(gs[3*i+1, :], sharex=msa_ax, aspect=block_rows*height_ratio/(data_max - data_min))

        block = im[:, i*sym_length*block_cols:(i+1)*sym_length*block_cols]
        x_left, x_right = x_start + i * block_cols, x_start + i * block_cols + block.shape[1] // sym_length
        msa_ax.imshow(block, extent=[x_left, x_right, block_rows, 0])
        msa_ax.set_yticks([y+0.5 for y, msa_label in zip(range(ROWS), msa_labels) if msa_label != ''])
        msa_ax.set_yticklabels([msa_label for msa_label in msa_labels if msa_label != ''])
        msa_ax.tick_params(axis='y', length=msa_ticklength, width=msa_tickwidth, pad=msa_tickpad, labelsize=msa_labelsize)
        msa_ax.xaxis.set_visible(False)
        for spine in ['left', 'right', 'top', 'bottom']:
            msa_ax.spines[spine].set_visible(False)

        for d in data:
            data_ax.plot(range(x_left, x_right), d[i*block_cols:i*block_cols + block.shape[1]//sym_length])
        data_ax.set_ylim(data_min, data_max)
        data_ax.tick_params(axis='y', labelsize=y_labelsize)
        data_ax.tick_params(axis='x', labelsize=x_labelsize)

    # Draw legend
    handles = []
    if msa_legend:
        color2syms = {}
        for sym, color in sym2color.items():
            if sym in gap2color:
                continue
            try:
                color2syms[color].append(sym)
            except KeyError:
                color2syms[color] = [sym]
        for color, syms in color2syms.items():
            handles.append(Line2D([], [], color=f'#{color}', linestyle='None', marker='.', label=', '.join(syms)))
    if data_labels:
        if handles:  # Add spacer entry if MSA legend is also drawn
            handles.append(Line2D([], [], color='none'))
        colors = rcParams['axes.prop_cycle'].by_key()['color']
        for i, data_label in enumerate(data_labels):
            handles.append(Line2D([], [], color=colors[i % len(colors)], label=data_label))
    if handles:
        fig.legend(handles=handles, **legend_kwargs)
    return fig


def plot_msa(msa, figsize=(12, 6),
             msa_labels=None, msa_labelsize=6, msa_length=0, msa_width=0.5, msa_pad=1,
             x_start=0, x_labelsize=6,
             hspace=0.5, sym_length=7, sym_height=7,
             msa_legend=False, legend_kwargs=None,
             block_cols=None, sym2color=None, gap2color=None):
    """Plot MSA as matplotlib figure.

    Parameters
    ----------
    msa: list of strings
    figsize: 2-tuple of floats
        Figure dimension (width, height) in inches.
    msa_labels: list of strings
        Strings of labels of sequences in MSA in the same order as in msa.
        len(msa) must match len(msa_labels). Empty strings are not displayed.
    msa_labelsize: float
        Font size of MSA labels.
    msa_length: float
        Length of MSA ticks.
    msa_width: float
        Width of MSA ticks.
    msa_pad: float
        Padding between MSA ticks and labels.
    x_start: int
        Starting value of x-axis.
    x_labelsize: float
        Font size of x-axis labels.
    hspace: float
        Padding between blocks of MSA as fraction of height of block.
    sym_length: int
        Number of pixels in length of the rectangles for each symbol.
    sym_height: int
        Number of pixels in height of the rectangles for each symbol.
    msa_legend: bool
        True if legend is drawn.
    legend_kwargs: dict
        Additional kwargs passed to fig.legend call.
    block_cols: int
        The number of columns in each block. Overrides determination of optimal
        number if not None.
    sym2color: dict
        Mapping of symbols to colors given as hex strings. These values
        determine the background color of the rectangles.
    gap2color: dict
        Mapping of symbols to colors given as hex strings. These values
        determine the foreground color of the dash in gap symbols.

    Returns
    -------
    fig: Figure (matplotlib)
    """
    # Define functions and globals
    ROWS, COLS = len(msa), len(msa[0])
    ASPECT = figsize[0] / figsize[1]

    def get_dims(block_cols, ROWS, COLS):
        plot_length = block_cols  # Length of final plot
        block_num = COLS // block_cols - (1 if COLS % block_cols == 0 else 0)  # Number of blocks in addition to the first
        plot_height = (1 + hspace) * ROWS * block_num + ROWS  # Height of final image
        return plot_length, plot_height

    # Set options
    if msa_labels is None:
        msa_labels = ROWS * ['']
    if len(msa_labels) != ROWS:
        raise RuntimeError('len(msa_labels) does not match len(msa)')
    if legend_kwargs is None:
        legend_kwargs = {}
    if block_cols is None:
        block_cols = _get_block_cols(ROWS, COLS, ASPECT, get_dims)
    if sym2color is None:
        sym2color = default_sym2color
    if gap2color is None:
        gap2color = default_gap2color
    block_num = COLS // block_cols + (1 if COLS % block_cols > 0 else 0)  # Number of blocks
    block_rows = len(msa)

    # Draw axes
    im = draw_msa(msa, block_cols=len(msa[0]), sym_length=sym_length, sym_height=sym_height, sym2color=sym2color, gap2color=gap2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(block_num, 1, figure=fig, hspace=hspace)
    for i in range(block_num):
        msa_ax = fig.add_subplot(gs[i, :])

        block = im[:, i*sym_length*block_cols:(i+1)*sym_length*block_cols]
        x_left, x_right = x_start + i * block_cols, x_start + i * block_cols + block.shape[1] // sym_length
        msa_ax.imshow(block, extent=[x_left, x_right, block_rows, 0])
        msa_ax.set_yticks([y+0.5 for y, msa_label in zip(range(ROWS), msa_labels) if msa_label != ''])
        msa_ax.set_yticklabels([msa_label for msa_label in msa_labels if msa_label != ''])
        msa_ax.tick_params(axis='y', length=msa_length, width=msa_width, pad=msa_pad, labelsize=msa_labelsize)
        msa_ax.tick_params(axis='x', labelsize=x_labelsize)
        for spine in ['left', 'right', 'top', 'bottom']:
            msa_ax.spines[spine].set_visible(False)

    # Draw legend
    handles = []
    if msa_legend:
        color2syms = {}
        for sym, color in sym2color.items():
            if sym in gap2color:
                continue
            try:
                color2syms[color].append(sym)
            except KeyError:
                color2syms[color] = [sym]
        for color, syms in color2syms.items():
            handles.append(Line2D([], [], color=f'#{color}', linestyle='None', marker='.', label=', '.join(syms)))
    if handles:
        fig.legend(handles=handles, **legend_kwargs)
    return fig


def plot_tree(tree, tip_labels=True, support_labels=False,
              color=None, linewidth=None,
              tip_fontsize=None, tip_fontcolor=None, tip_offset=0,
              support_fontsize=None, support_fontcolor=None,
              support_ha='center', support_va='top',
              support_hoffset=0, support_voffset=0):
    """Draw tree using matplotlib.

    Parameters
    ----------
    tree: TreeNode (skbio)
    tip_labels: bool
        Toggle drawing tip labels.
    support_labels: bool
        Toggle drawing support labels. Support labels are obtained from
        support attribute. Use assign_supports to extract the numerical
        values from the node labels.
    color: color (matplotlib)
    linewidth: float
    tip_fontsize: float or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
    tip_fontcolor: color (matplotlib)
    tip_offset: float
        The offset of the tip label from the end of the branch tip.
    support_fontsize: float
    support_fontcolor: color (matplotlib)
    support_ha: {'left', 'right', 'center'}
        The horizontal alignment of the support label relative to the branch.
    support_va: {'center', 'top', 'bottom', 'baseline', 'center_baseline'}
        The vertical alignment of the support label relative to the branch.
    support_voffset: float
        The vertical offset of the support label relative to its alignment.
    support_hoffset: float
        The horizontal offset of the support label relative to its alignment.

    Returns
    -------
    fig: Figure (matplotlib)
    ax: Axes (matplotlib)
    """
    # Set options
    if color is None:
        color = 'black'
    if linewidth is None:
        linewidth = rcParams['lines.linewidth']
    if tip_fontsize is None:
        tip_fontsize = rcParams['font.size']
    if tip_fontcolor is None:
        tip_fontcolor = 'black'
    if support_fontsize is None:
        support_fontsize = rcParams['font.size']
    if support_fontcolor is None:
        support_fontcolor = 'black'
    if support_ha == 'left':
        get_x = lambda node: x_pos.get(node.parent, 0)  # In case root node
    elif support_ha == 'center':
        get_x = lambda node: (x_pos.get(node.parent, 0) + x_pos[node]) / 2
    elif support_ha == 'right':
        get_x = lambda node: x_pos[node]
    else:
        raise ValueError(f"'{support_ha}' is not a valid value for support_ha; "
                         f"supported values are 'left', 'right', 'center'")

    # Make mapping of each node to its horizontal positions
    x_pos = {}
    for node in tree.preorder():  # Return nodes on the way in
        length = node.length if node.length is not None else 0
        depth = x_pos[node.parent] if node.parent else 0  # Checks for root node
        x_pos[node] = depth + length

    # Make mapping of each node to its vertical position
    tips = list(tree.tips())
    y_max = len(tips)
    y_pos = {tip: y_max - i for i, tip in enumerate(reversed(tips))}
    for node in tree.postorder():  # Return nodes on the way out
        if node.children:
            y_pos[node] = (y_pos[node.children[0]] + y_pos[node.children[-1]]) / 2

    # Plot lines and text
    lines = []
    fig, ax = plt.subplots()
    for node in tree.traverse():
        if node.parent:  # Horizontal line of node
            x0, x1 = x_pos[node.parent], x_pos[node]
            y = y_pos[node]
            lines.append(([(x0, y), (x1, y)], linewidth, color))
        if node.children:  # Vertical line of node
            x = x_pos[node]
            y0, y1 = y_pos[node.children[-1]], y_pos[node.children[0]]
            lines.append(([(x, y0), (x, y1)], linewidth, color))
        if tip_labels and node.is_tip():  # Write tip names
            ax.text(x_pos[node] + tip_offset, y_pos[node], node.name, verticalalignment='center',
                    fontsize=tip_fontsize, color=tip_fontcolor)
        if support_labels and not node.is_tip():  # Write support values
            ax.text(get_x(node) + support_hoffset, y_pos[node] + support_voffset, node.support,
                    fontsize=support_fontsize, color=support_fontcolor,
                    horizontalalignment=support_ha, verticalalignment=support_va)
    lc_args = {key: value for key, value in zip(['segments', 'linewidth', 'color'], zip(*lines))}
    ax.add_collection(LineCollection(**lc_args))

    # Adjust axes and add labels
    x_max = max(x_pos.values())
    y_max = max(y_pos.values())
    ax.set_xlim(-0.05 * x_max, 1.25 * x_max)  # Add margins around the tree to prevent overlapping the axes
    ax.set_ylim(y_max + 0.8, 0.2)  # Also invert the y-axis (origin at the top)
    ax.set_xlabel('branch length')

    return fig, ax


default_sym2color = {'A': '6DD7A1', 'I': '55C08C', 'L': '55C08C', 'V': '55C08C', 'M': '55C08C',
                     'F': 'B897EC', 'Y': 'B897EC', 'W': 'A180D2',
                     'S': 'FFBE74', 'T': 'FFBE74',
                     'N': '77EAF4', 'Q': '77EAF4',
                     'D': 'EE8485', 'E': 'EE8485',
                     'H': '96C4FF', 'K': '7FADEA', 'R': '7FADEA',
                     'C': 'FAED70', 'G': 'E2DEDD', 'P': 'FFB1F1',
                     'X': '93908F', '-': 'FFFFFF', '.': '3F3F3F'}
default_gap2color = {'-': '3F3F3F', '.': 'FFFFFF'}
