"""Drawing functions for specialized data."""

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import rcParams


def draw_tree(tree, tip_labels=True, support_labels=False,
         color=None, linewidth=None,
         tip_fontsize=None, tip_fontcolor=None, tip_offset=0,
         support_fontsize=None, support_fontcolor=None,
         support_ha='center', support_va='top',
         support_hoffset=0, support_voffset=0):
    """Draw tree using matplotlib.

    Parameters
    ----------
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
