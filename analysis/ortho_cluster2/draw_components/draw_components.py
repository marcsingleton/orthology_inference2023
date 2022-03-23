"""Draw largest connected components."""

import os
from math import exp

import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace

# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':') for adj in adjs.split(',')]

# Load connected components
components = {}
with open('../connect_hit_graph/out/components.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, nodes = line.rstrip().split('\t')
        components[component_id] = set(nodes.split(','))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

component_records = sorted(components.items(), key=lambda x: len(x[1]), reverse=True)
for i, (component_id, component) in enumerate(component_records[:50]):  # 50 largest components
    subgraph = {node: graph[node] for node in component}

    # Create graph
    nx_graph = nx.Graph()
    for node, adjs in subgraph.items():
        nx_graph.add_node(node)
        for adj, w in adjs:
            if (node, adj) in nx_graph.edges:
                edge_data = nx_graph.get_edge_data(node, adj)
                edge_data['weight'] = edge_data.get('weight', 0) + float(w)  # Sum edge weights
            else:
                nx_graph.add_edge(node, adj, weight=float(w))

    # Get positions and canvas limits
    positions = nx.kamada_kawai_layout(nx_graph, weight=None)  # Weight is None as otherwise is used for layout
    xs = [xy[0] for xy in positions.values()]
    xmin, xmax = min(xs), max(xs)
    xlen = xmax - xmin
    ys = [xy[1] for xy in positions.values()]
    ymin, ymax = min(ys), max(ys)
    ylen = ymax - ymin

    # Adjust dimensions so aspect ratio is 1:1
    max_dim = 8
    if xlen > ylen:
        figsize = (max_dim, max_dim * ylen/xlen)
    else:
        figsize = (max_dim * xlen/ylen, max_dim)

    # Determine best position for legend
    locs = ['lower left', 'lower right', 'upper left', 'upper right']
    for x, y in positions.values():
        xnorm = (x-xmin)/(xmax-xmin)
        ynorm = (y-ymin)/(ymax-ymin)
        if xnorm > 0.85 and ynorm > 0.95:
            locs.remove('upper right')
        elif xnorm < 0.15 and ynorm > 0.95:
            locs.remove('upper left')
        elif xnorm > 0.85 and ynorm < 0.05:
            locs.remove('lower right')
        elif xnorm < 0.15 and ynorm < 0.05:
            locs.remove('lower left')

    # Draw graph labeled by source
    node_size = 25 / (1 + exp(0.01 * (len(subgraph) - 400))) + 10  # Adjust node size
    FB = [node for node in nx_graph.nodes if node.startswith('FBpp')]
    NCBI = nx_graph.nodes - FB

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    nx.draw_networkx_edges(nx_graph, positions, alpha=0.25, width=0.5)
    nx.draw_networkx_nodes(NCBI, positions, node_size=node_size, linewidths=0, node_color='C0', label='NCBI')
    nx.draw_networkx_nodes(FB, positions, node_size=node_size, linewidths=0, node_color='C1', label='FlyBase')

    fig.legend(markerscale=(1 if node_size > 22.5 else 22.5/node_size), loc=locs[-1])
    fig.tight_layout()
    ax.axis('off')
    fig.savefig(f'out/{i:02}-{component_id}_source.png')
    plt.close()

    # Draw graph labeled by edge
    node_size = 20 / (1 + exp(0.01 * (len(subgraph) - 400))) + 10  # Adjust node size
    edges = sorted(nx_graph.edges, key=lambda x: nx_graph.get_edge_data(*x)['weight'])

    ws = [nx_graph.get_edge_data(*edge)['weight'] for edge in edges]
    cmap0 = mpl.cm.get_cmap('magma_r', 320)
    cmap1 = mpl.colors.ListedColormap(cmap0(linspace(0.25, 1, 256)))
    norm = mpl.colors.Normalize(min(ws), max(ws))

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    nx.draw_networkx_edges(nx_graph, positions, edgelist=edges, alpha=0.375, width=0.75, edge_color=ws, edge_cmap=cmap1)
    nx.draw_networkx_nodes(nx_graph, positions, node_size=node_size, linewidths=0, node_color='#333333')

    cax = inset_axes(ax, width=1, height=0.1, loc=locs[-1])
    dws = max(ws) - min(ws)
    ticks = [10*round((min(ws)+dws/4) / 10), 10*round((max(ws)-dws/4) / 10)]
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap1), cax=cax, ticks=ticks, orientation='horizontal')

    fig.tight_layout()
    ax.axis('off')
    fig.savefig(f'out/{i:02}-{component_id}_edge.png')
    plt.close()

"""
DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""