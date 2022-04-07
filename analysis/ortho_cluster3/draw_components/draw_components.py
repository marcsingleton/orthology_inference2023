"""Draw largest connected components."""

import os
from math import exp

import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace


def load_OGs(path):
    OGs = {}
    with open(path) as file:
        file.readline()  # Skip header
        for line in file:
            component_id, _, _, edges = line.rstrip('\n').split('\t')
            gnids = {node for edge in edges.split(',') for node in edge.split(':')}
            try:
                OGs[component_id].append(gnids)
            except KeyError:
                OGs[component_id] = [gnids]
    return OGs


def blend_colors(colors):
    rgbs = [[int(color[i:i+2], 16) for i in range(0, 6, 2)] for color in colors]
    avg = [int(sum(c) / len(c)) for c in zip(*rgbs)]
    return '#' + ''.join([hex(c)[2:] for c in avg])


def get_node_colors(graph, OGs):
    cycle = ['4E79A7', 'F28E2B', 'E15759', '76B7B2', '59A14F', 'EDC948', 'B07AA1', 'FF9DA7', '9C755F', 'BAB0AC']
    node2colors = {node: [] for node in graph.nodes}
    for i, OG in enumerate(OGs):
        for node in OG:
            node2colors[node].append(cycle[i % 10])

    node_colors = []
    for node in graph.nodes:
        colors = node2colors[node]
        if colors:
            node_colors.append(blend_colors(colors))
        else:
            node_colors.append('#1B1B1B')
    return node_colors


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
        component_id, nodes = line.rstrip('\n').split('\t')
        components[component_id] = set(nodes.split(','))

# Load OGs
OG3s = load_OGs('../cluster3_graph/out/clusters.tsv')
OG4s = load_OGs('../cluster4+_graph/out/4clique/clusters.tsv')
OG5s = load_OGs('../cluster4+_graph/out/5clique/clusters.tsv')
OG6s = load_OGs('../cluster4+_graph/out/6clique/clusters.tsv')

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
        if xnorm > 0.85 and ynorm > 0.95 and 'upper right' in locs:
            locs.remove('upper right')
        elif xnorm < 0.15 and ynorm > 0.95 and 'upper left' in locs:
            locs.remove('upper left')
        elif xnorm > 0.85 and ynorm < 0.05 and 'lower right' in locs:
            locs.remove('lower right')
        elif xnorm < 0.15 and ynorm < 0.05 and 'lower left' in locs:
            locs.remove('lower left')

    # Draw graph labeled by source
    node_size = 25 / (1 + exp(0.01 * (len(subgraph) - 400))) + 10  # Adjust node size
    FB = {node for node in nx_graph.nodes if node.startswith('FBpp')}
    NCBI = nx_graph.nodes - FB

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes((0.01, 0.01, 0.98, 0.98))
    nx.draw_networkx_edges(nx_graph, positions, alpha=0.25, width=0.5)
    nx.draw_networkx_nodes(NCBI, positions, node_size=node_size, linewidths=0, node_color='C0', label='NCBI')
    nx.draw_networkx_nodes(FB, positions, node_size=node_size, linewidths=0, node_color='C1', label='FlyBase')

    fig.legend(markerscale=(1 if node_size > 22.5 else 22.5/node_size), loc=locs[-1])
    ax.axis('off')
    fig.savefig(f'out/{i:02}_{component_id}_source.png', dpi=400)
    plt.close()

    # Draw graph labeled by cluster
    for k, OGks in zip(range(3, 7), [OG3s, OG4s, OG5s, OG6s]):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes((0.01, 0.01, 0.98, 0.98))
        nx.draw_networkx_edges(nx_graph, positions, alpha=0.25, width=0.5)
        nx.draw_networkx_nodes(nx_graph.nodes, positions, node_size=node_size, linewidths=0, node_color=get_node_colors(nx_graph, OGks.get(component_id, [])))

        ax.axis('off')
        fig.savefig(f'out/{i:02}_{component_id}_OG{k}.png', dpi=400)
        plt.close()

    # Draw graph labeled by edge
    node_size = 20 / (1 + exp(0.01 * (len(subgraph) - 400))) + 10  # Adjust node size
    edges = sorted(nx_graph.edges, key=lambda x: nx_graph.get_edge_data(*x)['weight'])

    ws = [nx_graph.get_edge_data(*edge)['weight'] for edge in edges]
    cmap0 = mpl.cm.get_cmap('magma_r', 320)
    cmap1 = mpl.colors.ListedColormap(cmap0(linspace(0.25, 1, 256)))
    norm = mpl.colors.Normalize(min(ws), max(ws))

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes((0.01, 0.01, 0.98, 0.98))
    nx.draw_networkx_edges(nx_graph, positions, edgelist=edges, alpha=0.375, width=0.75, edge_color=ws, edge_cmap=cmap1)
    nx.draw_networkx_nodes(nx_graph, positions, node_size=node_size, linewidths=0, node_color='#333333')

    cax = inset_axes(ax, width=1, height=0.1, loc=locs[-1])
    dws = max(ws) - min(ws)
    ticks = [10*round((min(ws)+dws/4) / 10), 10*round((max(ws)-dws/4) / 10)]
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap1), cax=cax, ticks=ticks, orientation='horizontal')

    ax.axis('off')
    fig.savefig(f'out/{i:02}_{component_id}_edge.png', dpi=400)
    plt.close()

"""
DEPENDENCIES
../cluster3_graph/cluster3_graph.py
    ../cluster3_graph/out/clusters.tsv
../cluster4+_graph/cluster4+_graph.py
    ../cluster4+_graph/out/4clique/clusters.tsv
    ../cluster4+_graph/out/5clique/clusters.tsv
    ../cluster4+_graph/out/6clique/clusters.tsv
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""