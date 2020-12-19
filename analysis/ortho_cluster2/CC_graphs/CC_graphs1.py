"""Draw largest connected components."""

import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import os
from math import exp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace

# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph1.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = [adj.split(':') for adj in adjs.split(',')]

# Load connected components
CCs = {}
with open('../connect_ggraph/out/gconnect1.txt') as file:
    for line in file:
        CCid, nodes = line.rstrip().split(':')
        CCs[CCid] = set(nodes.split(','))

# Make output directory
if not os.path.exists('out/ggraph1/'):
    os.makedirs('out/ggraph1/')  # Recursive folder creation

CCids = sorted(CCs, key=lambda x: len(CCs[x]), reverse=True)
for i, CCid in enumerate(CCids[:50]):  # 50 largest CCs
    subggraph = {node: ggraph[node] for node in CCs[CCid]}

    # Create graph and segment nodes by data source
    G = nx.Graph()
    for node, adjs in subggraph.items():
        G.add_node(node)
        for adj, w in adjs:
            if (node, adj) in G.edges:
                edge_data = G.get_edge_data(node, adj)
                edge_data['weight'] = edge_data.get('weight', 0) + int(w)  # Sum edge weights
            else:
                G.add_edge(node, adj, weight=int(w))

    # Get positions and canvas limits
    pos = nx.kamada_kawai_layout(G, weight=None)  # Weight is None as otherwise is used for layout
    xs = [xy[0] for xy in pos.values()]
    xmin, xmax = min(xs), max(xs)
    xlen = xmax - xmin
    ys = [xy[1] for xy in pos.values()]
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
    delset = set()
    for x, y in pos.values():
        xnorm = (x-xmin)/(xmax-xmin)
        ynorm = (y-ymin)/(ymax-ymin)
        if xnorm > 0.85 and ynorm > 0.95:
            delset.add('upper right')
        elif xnorm < 0.15 and ynorm > 0.95:
            delset.add('upper left')
        elif xnorm > 0.85 and ynorm < 0.05:
            delset.add('lower right')
        elif xnorm < 0.15 and ynorm < 0.05:
            delset.add('lower left')
    for loc in delset:
        locs.remove(loc)

    # Draw graph labeled by source
    node_size = 25/(1 + exp(0.01*(len(subggraph)-400))) + 10  # Adjust node size
    FB = [node for node in G.nodes if node.startswith('FBgn')]
    NCBI = G.nodes - FB

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.5)
    nx.draw_networkx_nodes(NCBI, pos, node_size=node_size, linewidths=0, node_color='C0', label='NCBI')
    nx.draw_networkx_nodes(FB, pos, node_size=node_size, linewidths=0, node_color='C1', label='FlyBase')

    fig.legend(markerscale=(1 if node_size > 22.5 else 22.5/node_size), loc=locs[-1])
    fig.tight_layout()
    ax.axis('off')
    fig.savefig(f'out/ggraph1/{i}_{CCid}_source.png')
    plt.close()

    # Draw graph labeled by edge
    node_size = 20/(1 + exp(0.01*(len(subggraph)-400))) + 10  # Adjust node size
    edges = sorted(G.edges, key=lambda x: G.get_edge_data(*x)['weight'])

    ws = [G.get_edge_data(*edge)['weight'] for edge in edges]
    cmap0 = mpl.cm.get_cmap('magma_r', 320)
    cmap1 = mpl.colors.ListedColormap(cmap0(linspace(0.25, 1, 256)))
    norm = mpl.colors.Normalize(min(ws), max(ws))

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    nx.draw_networkx_edges(G, pos, edgelist=edges, alpha=0.375, width=0.75, edge_color=ws, edge_cmap=cmap1)
    nx.draw_networkx_nodes(G, pos, node_size=node_size, linewidths=0, node_color='#333333')

    cax = inset_axes(ax, width=1, height=0.1, loc=locs[-1])
    dws = max(ws) - min(ws)
    ticks = [10*round((min(ws)+dws/4) / 10), 10*round((max(ws)-dws/4) / 10)]
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap1), cax=cax, ticks=ticks, orientation='horizontal')

    fig.tight_layout()
    ax.axis('off')
    fig.savefig(f'out/ggraph1/{i}_{CCid}_edge.png')
    plt.close()

"""
DEPENDENCIES
../connect_ggraph/connect_ggraph1.py
    ../connect_ggraph/out/gconnect1.txt
../hits2ggraph/hits2ggraph1.py
    ../hits2ggraph/out/ggraph1.tsv
"""