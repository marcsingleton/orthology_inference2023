"""Drawing functions for specialized data."""

from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory


def _get_aspect(block_columns, get_dimensions):
    length, height = get_dimensions(block_columns)
    return length / height


def _get_block_columns(msa_columns, aspect, get_dimensions):
    # Use binary search to find interval containing optimal block_columns
    # If an interval endpoint exactly equals the ratio, choose that endpoint
    # Otherwise choose the interval where the deviation changes sign between the two endpoints
    interval = (1, msa_columns)
    while interval[1] - interval[0] > 1:
        i0, i2 = interval
        i1 = (i0 + i2) // 2
        deltas = [_get_aspect(i, get_dimensions) - aspect for i in (i0, i1, i2)]

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
    block_columns = min(range(interval[0], interval[1] + 1), key=lambda x: abs(_get_aspect(x, get_dimensions) - aspect))  # Choose value that minimizes difference

    # Ensure last block is at least 50% of block_columns
    if msa_columns % block_columns < 0.5 * block_columns:  # Guarantees at least two blocks
        blocks_im = msa_columns // block_columns  # Total blocks minus 1
        block_columns += ceil((msa_columns % block_columns) / blocks_im)  # Distribute excess to other blocks
    return block_columns


def draw_msa(msa, *,
             aspect=2.5, hspace=25, sym_length=7, sym_height=7,
             block_columns=None, sym2color=None, gap2color=None):
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
    block_columns: int
        Number of columns in each block. Will override aspect if is not None.
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
    msa_rows, msa_columns = len(msa), len(msa[0])

    def get_dimensions(block_columns):
        im_length = sym_length * block_columns  # Length of final image
        block_number = ceil(msa_columns / block_columns)
        im_height = sym_height * msa_rows + (block_number - 1) * (sym_height * msa_rows + hspace)
        return im_length, im_height

    # Set options
    if block_columns is None:
        block_columns = _get_block_columns(msa_columns, aspect, get_dimensions)
    if sym2color is None:
        sym2color = default_sym2color
    if gap2color is None:
        gap2color = default_gap2color

    # Instantiate array and fill with values
    im_length, im_height = get_dimensions(block_columns)
    im = np.full((im_height, im_length, 3), 255, dtype='uint8')
    for i, seq in enumerate(msa):
        for j, sym in enumerate(seq):
            # Position of symbol rectangle
            block = j // block_columns
            y = (sym_height * msa_rows + hspace) * block + sym_height * i
            x = j % block_columns * sym_length

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


def plot_msa_data(msa, data, *,
                  msa_labels=None, msa_labelsize=6, msa_ticklength=0, msa_tickwidth=0.5, msa_tickpad=1,
                  x_start=0, x_labelsize=6, y_labelsize=6,
                  tree=None, tree_position=0, tree_width=0.1, tree_kwargs=None,
                  height_ratio=1, hspace=0.25, sym_length=7, sym_height=7,
                  figsize=(15, 6), left=0.05, right=0.95, top=0.95, bottom=0.05, anchor=(0.5, 0.5),
                  data_min=None, data_max=None,
                  data_labels=None, data_colors=None, data_alphas=None,
                  msa_legend=False, legend_kwargs=None,
                  block_columns=None, sym2color=None, gap2color=None):
    """Plot MSA with associated positional data as matplotlib figure.

    Parameters
    ----------
    msa: list of strings
    data: list of lists or two-dimensional ndarray
        Data series to be plotted where outer list is the series and the inner
        list is a value associated with each column of the msa.
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
    tree: TreeNode (skbio)
    tree_position: float
        Position of left edge of tree in figure coordinates.
    tree_width: float
        Width of tree in figure coordinates
    tree_kwargs: dict
        Dictionary of additional keyword arguments passed to plot_tree.
    height_ratio: float
        Height of data axes as fraction of height of block.
    hspace: float
        Padding between each data axes and subsequent block of MSA as fraction
        of average height of blocks and data axes.
    sym_length: int
        Number of pixels in length of the rectangles for each symbol.
    sym_height: int
        Number of pixels in height of the rectangles for each symbol.
    figsize: 2-tuple of floats
        Figure dimension (width, height) in inches.
    left, right, top, bottom: float
        Extent of the subplots as a fraction of figure width or height.
    anchor: (float, float)
        Anchor for subplots.
    data_min: float
        Minimum of y-axis across all data axes.
    data_max: float
        Maximum of y-axis across all data axes.
    data_labels: list of strings
        Labels for data series. If not None, will draw legend. Length must
        match the number of series in data.
    data_colors: list of colors
        Colors for data series. Length must match number of series in data.
    data_alphas: list of float
        Alphas for data series. Length must match number of series in data.
    msa_legend: bool
        True if MSA legend is drawn.
    legend_kwargs: dict
        Additional kwargs passed to fig.legend call.
    block_columns: int
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
    msa_rows, msa_columns = len(msa), len(msa[0])
    aspect = (figsize[0] * (right - left)) / (figsize[1] * (top - bottom))

    def get_dimensions(block_columns):
        plot_length = block_columns
        block_number = ceil(msa_columns / block_columns)
        hspace_effective = hspace * (1 + height_ratio) / 2
        plot_height = (1 + hspace_effective + height_ratio) * msa_rows + (block_number - 1) * (1 + 2 * hspace_effective + height_ratio) * msa_rows
        return plot_length, plot_height

    # Set options
    if msa_labels is None:
        msa_labels = msa_rows * ['']
    elif len(msa_labels) != msa_rows:
        raise ValueError('len(msa_labels) does not match len(msa)')
    if tree_kwargs is None:
        tree_kwargs = {}
    if legend_kwargs is None:
        legend_kwargs = {}
    if block_columns is None:
        block_columns = _get_block_columns(msa_columns, aspect, get_dimensions)
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
    if data_labels is not None and len(data_labels) != len(data):
        raise ValueError('len(data_labels) does not match len(data)')
    if data_colors is None:
        color_cycle = rcParams['axes.prop_cycle'].by_key()['color']
        data_colors = [color_cycle[i % len(color_cycle)] for i in range(len(data))]
    elif len(data_colors) != len(data):
        raise ValueError('len(data_colors) does not match len(data)')
    if data_alphas is None:
        data_alphas = [1 for _ in range(len(data))]
    elif len(data_alphas) != len(data):
        raise ValueError('len(data_alphas) does not match len(data)')
    block_number = ceil(msa_columns / block_columns)
    block_rows = len(msa)

    # Draw axes
    height_ratios = []
    for i in range(2*block_number):
        r = i % 2
        if r == 0:
            h = 1
        else:
            h = height_ratio
        height_ratios.append(h)

    im = draw_msa(msa, block_columns=len(msa[0]), sym_length=sym_length, sym_height=sym_height, sym2color=sym2color, gap2color=gap2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2*block_number, 1, figure=fig, left=left, right=right, top=top, bottom=bottom,
                  hspace=hspace, height_ratios=height_ratios)
    for i in range(block_number):
        msa_ax = fig.add_subplot(gs[2*i, :], anchor=anchor)
        data_ax = fig.add_subplot(gs[2*i+1, :], sharex=msa_ax, aspect=block_rows*height_ratio/(data_max - data_min), anchor=anchor)

        block = im[:, i * sym_length * block_columns:(i + 1) * sym_length * block_columns]
        x_left, x_right = x_start + i * block_columns, x_start + i * block_columns + block.shape[1] // sym_length
        msa_ax.imshow(block, extent=[x_left, x_right, block_rows, 0])
        msa_ax.set_yticks([y+0.5 for y, msa_label in zip(range(msa_rows), msa_labels) if msa_label != ''])
        msa_ax.set_yticklabels([msa_label for msa_label in msa_labels if msa_label != ''])
        msa_ax.tick_params(axis='y', length=msa_ticklength, width=msa_tickwidth, pad=msa_tickpad, labelsize=msa_labelsize)
        msa_ax.xaxis.set_visible(False)
        for spine in ['left', 'right', 'top', 'bottom']:
            msa_ax.spines[spine].set_visible(False)

        for series, color, alpha in zip(data, data_colors, data_alphas):
            data_ax.plot(range(x_left, x_right), series[i * block_columns:i * block_columns + block.shape[1] // sym_length],
                         color=color, alpha=alpha)
        data_ax.set_ylim(data_min, data_max)
        data_ax.tick_params(axis='y', labelsize=y_labelsize)
        data_ax.tick_params(axis='x', labelsize=x_labelsize)

        if tree:
            transform = blended_transform_factory(fig.transFigure, msa_ax.transAxes)
            tree_ax = msa_ax.inset_axes([tree_position, 0, tree_width, 1], transform=transform)
            plot_tree(tree, ax=tree_ax, ymin_pad=1 / (2 * msa_rows), ymax_pad=1 / (2 * msa_rows), **tree_kwargs)
            tree_ax.axis('off')

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
        for data_label, data_color in zip(data_labels, data_colors):
            handles.append(Line2D([], [], color=data_color, label=data_label))
    if handles:
        fig.legend(handles=handles, **legend_kwargs)
    return fig


def plot_msa(msa, *,
             msa_labels=None, msa_labelsize=6, msa_length=0, msa_width=0.5, msa_pad=1,
             x_start=0, x_labelsize=6,
             tree=None, tree_position=0, tree_width=0.1, tree_kwargs=None,
             hspace=0.25, sym_length=7, sym_height=7,
             figsize=(15, 6), left=0.05, right=0.95, top=0.95, bottom=0.05, anchor=(0.5, 0.5),
             msa_legend=False, legend_kwargs=None,
             block_columns=None, sym2color=None, gap2color=None):
    """Plot MSA as matplotlib figure.

    Parameters
    ----------
    msa: list of strings
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
    tree: TreeNode (skbio)
    tree_position: float
        Position of left edge of tree in figure coordinates.
    tree_width: float
        Width of tree in figure coordinates
    tree_kwargs: dict
        Dictionary of additional keyword arguments passed to plot_tree.
    hspace: float
        Padding between blocks of MSA as fraction of height of blocks.
    sym_length: int
        Number of pixels in length of the rectangles for each symbol.
    sym_height: int
        Number of pixels in height of the rectangles for each symbol.
    figsize: 2-tuple of floats
        Figure dimension (width, height) in inches.
    left, right, top, bottom: float
        Extent of the subplots as a fraction of figure width or height.
    anchor: (float, float)
        Anchor for subplots.
    msa_legend: bool
        True if legend is drawn.
    legend_kwargs: dict
        Additional kwargs passed to fig.legend call.
    block_columns: int
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
    msa_rows, msa_columns = len(msa), len(msa[0])
    aspect = (figsize[0] * (right - left)) / (figsize[1] * (top - bottom))

    def get_dimensions(block_columns):
        plot_length = block_columns
        block_number = ceil(msa_columns / block_columns)
        plot_height = msa_rows + (block_number - 1) * (1 + hspace) * msa_rows
        return plot_length, plot_height

    # Set options
    if msa_labels is None:
        msa_labels = msa_rows * ['']
    elif len(msa_labels) != msa_rows:
        raise ValueError('len(msa_labels) does not match len(msa)')
    if tree_kwargs is None:
        tree_kwargs = {}
    if legend_kwargs is None:
        legend_kwargs = {}
    if block_columns is None:
        block_columns = _get_block_columns(msa_columns, aspect, get_dimensions)
    if sym2color is None:
        sym2color = default_sym2color
    if gap2color is None:
        gap2color = default_gap2color
    block_number = ceil(msa_columns / block_columns)
    block_rows = len(msa)

    # Draw axes
    im = draw_msa(msa, block_columns=len(msa[0]), sym_length=sym_length, sym_height=sym_height, sym2color=sym2color, gap2color=gap2color)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(block_number, 1, figure=fig, left=left, right=right, top=top, bottom=bottom, hspace=hspace)
    for i in range(block_number):
        msa_ax = fig.add_subplot(gs[i, :], anchor=anchor)

        block = im[:, i * sym_length * block_columns:(i + 1) * sym_length * block_columns]
        x_left, x_right = x_start + i * block_columns, x_start + i * block_columns + block.shape[1] // sym_length
        msa_ax.imshow(block, extent=[x_left, x_right, block_rows, 0])
        msa_ax.set_yticks([y + 0.5 for y, msa_label in zip(range(msa_rows), msa_labels) if msa_label != ''])
        msa_ax.set_yticklabels([msa_label for msa_label in msa_labels if msa_label != ''])
        msa_ax.tick_params(axis='y', length=msa_length, width=msa_width, pad=msa_pad, labelsize=msa_labelsize)
        msa_ax.tick_params(axis='x', labelsize=x_labelsize)
        for spine in ['left', 'right', 'top', 'bottom']:
            msa_ax.spines[spine].set_visible(False)

        if tree:
            transform = blended_transform_factory(fig.transFigure, msa_ax.transAxes)
            tree_ax = msa_ax.inset_axes([tree_position, 0, tree_width, 1], transform=transform)
            plot_tree(tree, ax=tree_ax, ymin_pad=1 / (2 * msa_rows), ymax_pad=1 / (2 * msa_rows), **tree_kwargs)
            tree_ax.axis('off')

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


def plot_tree(tree, *, ax=None,
              linecolor=None, linewidth=None,
              tip_labels=True, tip_fontsize=None, tip_fontcolor=None, tip_offset=0.005,
              support_labels=False, support_format_spec=None, support_fontsize=None, support_fontcolor=None,
              support_ha='center', support_va='top', support_hoffset=0, support_voffset=0,
              xmin_pad=0.01, xmax_pad=0.1, ymin_pad=0.025, ymax_pad=0.025):
    """Draw tree using matplotlib.

    Parameters
    ----------
    tree: TreeNode (skbio)
    ax: Axes (matplotlib)
        Axes used to draw tree. If None, a new Figure and Axes are created.
    linecolor: color (matplotlib) or dict
        If is dict, linecolor is a mapping of nodes to colors.
    linewidth: float
        Width of branches.
    tip_labels: bool
        Toggle drawing tip labels.
    tip_fontsize: float or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
    tip_fontcolor: color (matplotlib)
    tip_offset: float
        The offset of the tip label from the end of the branch tip in axes
        coordinates.
    support_labels: bool
        Toggle drawing support labels. Support labels are obtained from
        support attribute. If the attribute is None, it is ignored. Use
        assign_supports to extract the numerical values from the node labels.
    support_format_spec: str
        Format specification for supports using the format specification mini-
        language.
    support_fontsize: float
    support_fontcolor: color (matplotlib)
    support_ha: {'left', 'right', 'center'}
        The horizontal alignment of the support label relative to the branch.
    support_va: {'center', 'top', 'bottom', 'baseline', 'center_baseline'}
        The vertical alignment of the support label relative to the branch.
    support_voffset: float
        The vertical offset of the support label relative to its alignment in
        axes coordinates.
    support_hoffset: float
        The horizontal offset of the support label relative to its alignment in
        axes coordinates.
    xmin_pad, xmax_pad, ymin_pad, ymax_pad: float
        Fraction of respective data ranges to pad on lower and upper bounds of
        respective axes.

    Returns
    -------
    fig: Figure (matplotlib)
    ax: Axes (matplotlib)
    """
    # Set options
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    if not isinstance(linecolor, dict):
        if linecolor is None:
            linecolor = 'black'
        node2color = defaultdict(lambda x=linecolor: x)  # Use default parameter to ensure lambda captures value
    else:
        node2color = linecolor
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
        get_x = lambda node: xpos.get(node.parent, 0)  # In case root node
    elif support_ha == 'center':
        get_x = lambda node: (xpos.get(node.parent, 0) + xpos[node]) / 2
    elif support_ha == 'right':
        get_x = lambda node: xpos[node]
    else:
        raise ValueError(f"'{support_ha}' is not a valid value for support_ha; "
                         f"supported values are 'left', 'right', 'center'")

    # Make mapping of each node to its horizontal positions
    xpos = {}
    for node in tree.preorder():  # Return nodes on the way in
        length = node.length if node.length is not None else 0
        depth = xpos[node.parent] if node.parent else 0  # Checks for root node
        xpos[node] = depth + length

    # Make mapping of each node to its vertical position
    tips = list(tree.tips())
    ymax = len(tips)
    ypos = {tip: ymax - i for i, tip in enumerate(reversed(tips))}
    for node in tree.postorder():  # Return nodes on the way out
        if node.children:
            ypos[node] = (ypos[node.children[0]] + ypos[node.children[-1]]) / 2

    # Adjust axes and add labels
    xmin, xmax = min(xpos.values()), max(xpos.values())
    xrange = xmax - xmin
    ax.set_xlim(xmin - xmin_pad * xrange, xmax + xmax_pad * xrange)

    ymin, ymax = min(ypos.values()), max(ypos.values())
    yrange = ymax - ymin
    ax.set_ylim(ymax + ymax_pad * yrange, ymin - ymin_pad * yrange)  # Invert the y-axis (origin at the top)

    # Plot lines and text
    lines = []
    transform = ax.transLimits.inverted()  # Axes to data coordinates
    for node in tree.traverse():
        linecolor = node2color[node]
        if node.parent:  # Horizontal line of node
            x0, x1 = xpos[node.parent], xpos[node]
            y = ypos[node]
            lines.append(([(x0, y), (x1, y)], linewidth, linecolor))
        if node.children:  # Vertical line of node
            x = xpos[node]
            y0, y1 = ypos[node.children[-1]], ypos[node.children[0]]
            lines.append(([(x, y0), (x, y1)], linewidth, linecolor))
        if tip_labels and node.is_tip():  # Write tip names
            dx, _ = transform.transform((tip_offset, 0)) - transform.transform((0, 0))
            ax.text(xpos[node] + dx, ypos[node], node.name, verticalalignment='center',
                    fontsize=tip_fontsize, color=tip_fontcolor)
        if support_labels and node.support is not None and not node.is_tip():  # Write support values if not None and not tip
            if support_format_spec is None:
                support_string = str(node.support)
            else:
                support_string = f'{node.support:{support_format_spec}}'
            dx, dy = transform.transform((support_hoffset, support_voffset)) - transform.transform((0, 0))
            ax.text(get_x(node) + dx, ypos[node] + dy, support_string,
                    fontsize=support_fontsize, color=support_fontcolor,
                    horizontalalignment=support_ha, verticalalignment=support_va)
    lc_args = {key: value for key, value in zip(['segments', 'linewidth', 'color'], zip(*lines))}
    ax.add_collection(LineCollection(**lc_args, capstyle='projecting'))

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
