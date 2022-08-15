"""Draw largest connected components."""

import os
from math import floor, log10

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from numpy import linspace

num_components = 100

fig_width = 7.5
dpi = 400

margin = 0
margin_data = 0.1
axes_rect = (margin, margin, 1 - margin, 1 - margin)

node_color2 = '#444444'
node_slope2 = 125
node_intercept2 = 0
edge_alpha2 = 0.3
edge_width2 = 0.75

cmap_base = mpl.cm.get_cmap('plasma_r', 320)  # 64 additional levels so truncated cmap has 256
cmap = mpl.colors.ListedColormap(cmap_base(linspace(0.25, 1, 256)))
cbar_fontsize = 'medium'
cbar_width = 1.5
cbar_height = 0.1
cbar_offset = 0.01

# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':') for adj in adjs.split(',')]

# Load connected components
components = {}
with open('../connect_hit_graph/out/components.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        components[fields['component_id']] = set(fields['ppids'].split(','))

if not os.path.exists('out/'):
    os.mkdir('out/')

component_records = sorted(components.items(), key=lambda x: len(x[1]), reverse=True)
for i, (component_id, component) in enumerate(component_records[:num_components]):  # Largest components only
    # Create graph
    subgraph = {node: graph[node] for node in component}

    nx_graph = nx.Graph()
    for node, adjs in sorted(subgraph.items()):
        nx_graph.add_node(node)
        for adj, w in adjs:
            if (node, adj) in nx_graph.edges:
                nx_graph.edges[node, adj]['weight'] += float(w)
            else:
                nx_graph.add_edge(node, adj, weight=float(w))

    # Get positions and axes limits
    positions = nx.kamada_kawai_layout(nx_graph, weight=None)  # Weight is None as otherwise is used for layout
    xs = [xy[0] for xy in positions.values()]
    xmin, xmax = min(xs), max(xs)
    xlen = xmax - xmin
    ys = [xy[1] for xy in positions.values()]
    ymin, ymax = min(ys), max(ys)
    ylen = ymax - ymin

    # Rotate positions
    if xlen / ylen < 1:  # Make width longer side
        xmin, xmax, ymin, ymax = ymin, ymax, -xmax, -xmin
        xlen, ylen = ylen, xlen
        positions = {node: (y, -x) for node, (x, y) in positions.items()}

    # Adjust dimensions so aspect ratio is 1:1
    fig_height = fig_width * ylen / xlen
    figsize = (fig_width, fig_height)
    node_size2 = node_slope2 * fig_width * fig_height / len(subgraph) + node_intercept2

    # Draw graph labeled by edge
    edges = sorted(nx_graph.edges, key=lambda x: nx_graph.edges[x]['weight'])
    ws = [w for _, _, w in nx_graph.edges.data('weight')]
    wmin, wmax = min(ws), max(ws)
    wlen = wmax - wmin

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axes_rect)

    nx.draw_networkx_edges(nx_graph, positions, ax=ax, edgelist=edges, alpha=edge_alpha2, width=edge_width2, edge_color=ws, edge_cmap=cmap)
    nx.draw_networkx_nodes(nx_graph, positions, ax=ax, node_size=node_size2, linewidths=0, node_color=node_color2)
    ax.set_xlim((xmin - margin_data * xlen, xmax + margin_data * xlen))  # Set manually because draw_networkx_edges hard codes the data limits with 5% padding
    ax.set_ylim((ymin - margin_data * ylen, ymax + margin_data * ylen))

    cax = fig.add_axes((1 - cbar_width / fig_width - cbar_offset, 1 - cbar_height / fig_height - cbar_offset,
                        cbar_width / fig_width, cbar_height / fig_height))
    ticks = [wmin + wlen / 4, wmax - wlen / 4]
    ticklabels = [f'{round(tick, -floor(log10(tick)) + 2):5g}' for tick in ticks]  # Round to third significant digit of difference between ticks
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(wmin, wmax), cmap=cmap), cax=cax, orientation='horizontal')
    cbar.ax.set_xticks(ticks, ticklabels, fontsize=cbar_fontsize)

    ax.axis('off')
    fig.savefig(f'out/{i:02}_{component_id}_edge.png', dpi=dpi)
    plt.close()

"""
DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""