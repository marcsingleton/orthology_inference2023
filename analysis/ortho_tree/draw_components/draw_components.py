"""Draw largest connected components."""

import os
from math import floor, log10

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from numpy import linspace


def blend_colors(colors):
    rgbs = [[int(color[i:i+2], 16) for i in range(0, 6, 2)] for color in colors]
    avg = [int(sum(c) / len(c)) for c in zip(*rgbs)]
    return '#' + ''.join([f'{c:X}' for c in avg])


def get_node_colors(graph, OGs):
    cycle = ['4E79A7', 'F28E2B', 'E15759', '76B7B2', '59A14F', 'EDC948', 'B07AA1', 'FF9DA7', '9C755F', 'BAB0AC']
    node2colors = {node: [] for node in graph.nodes}
    for i, OG in enumerate(OGs):
        for node in OG:
            node2colors[node].append(cycle[i % len(cycle)])

    node_colors = []
    for node in graph.nodes:
        colors = node2colors[node]
        if colors:
            node_colors.append(blend_colors(colors))
        else:
            node_colors.append('#1B1B1B')
    return node_colors


num_components = 100

fig_width = 7.5
margin = 0.025
axes_rect = (margin, margin, 1 - margin, 1 - margin)
dpi = 400

cbar_width = 1.5
cbar_height = 0.1
cbar_margin = 0.01

node_slope = 75
node_intercept = 3

edge_alpha1 = 0.15
edge_width1 = 0.5
edge_alpha2 = 0.375
edge_width2 = 0.75

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
        xlen, ylen = ylen, xlen
        positions = {node: (y, -x) for node, (x, y) in positions.items()}

    # Adjust dimensions so aspect ratio is 1:1
    fig_height = fig_width * ylen / xlen
    figsize = (fig_width, fig_height)
    node_size = node_slope * fig_width * fig_height / len(subgraph) + node_intercept  # Node size is inversely proportional to node density

    # Draw graph labeled by source
    FB = {node for node in nx_graph.nodes if node.startswith('FBpp')}
    NCBI = nx_graph.nodes - FB

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axes_rect)
    nx.draw_networkx_edges(nx_graph, positions, ax=ax, alpha=0.25, width=0.5)
    nx.draw_networkx_nodes(NCBI, positions, ax=ax, node_size=node_size, linewidths=0, node_color='C0', label='NCBI')
    nx.draw_networkx_nodes(FB, positions, ax=ax, node_size=node_size, linewidths=0, node_color='C1', label='FlyBase')

    fig.legend(markerscale=25/node_size, loc='upper right')  # Node size on legend is constant
    ax.axis('off')
    fig.savefig(f'out/{i:02}_{component_id}_source.png', dpi=dpi)
    plt.close()

    # Draw graph labeled by edge
    edges = sorted(nx_graph.edges, key=lambda x: nx_graph.edges[x]['weight'])
    ws = [w for _, _, w in nx_graph.edges.data('weight')]
    wmin, wmax = min(ws), max(ws)
    wlen = wmax - wmin

    cmap0 = mpl.cm.get_cmap('magma_r', 320)  # 64 additional levels so truncated cmap has 256
    cmap1 = mpl.colors.ListedColormap(cmap0(linspace(0.25, 1, 256)))
    norm = mpl.colors.Normalize(wmin, wmax)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axes_rect)
    nx.draw_networkx_edges(nx_graph, positions, ax=ax, edgelist=edges, alpha=edge_alpha2, width=edge_width2, edge_color=ws, edge_cmap=cmap1)
    nx.draw_networkx_nodes(nx_graph, positions, ax=ax, node_size=node_size, linewidths=0, node_color='#444444')

    cax = fig.add_axes((1 - cbar_width/fig_width - cbar_margin, 1 - cbar_height/fig_height - cbar_margin,
                        cbar_width/fig_width, cbar_height/fig_height))
    ndigits = -floor(log10(wlen / 2)) + 1  # Round to second significant digit of difference between ticks
    ticks = [round(wmin + wlen / 4, ndigits), round(wmax - wlen / 4, ndigits)]
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap1), cax=cax, ticks=ticks, orientation='horizontal')

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