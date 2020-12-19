"""Draw largest connected components."""

import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import os
from math import exp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace


def load_OGs(path):
    OGs = {}
    with open(path) as file:
        for line in file:
            CCid, OGid, edges = line.rstrip().split(':')
            gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
            try:
                OGs[CCid].append(gnids)
            except KeyError:
                OGs[CCid] = [gnids]
    return OGs


def blend_colors(colors):
    rgbs = [[int(color[i:i+2], 16) for i in range(0, 6, 2)] for color in colors]
    avg = [int(sum(c) / len(c)) for c in zip(*rgbs)]
    return '#' + ''.join([hex(c)[2:] for c in avg])


def get_node_colors(graph, OGs):
    cycle = ['4e79a7', 'f28e2b', 'e15759', '76b7b2', '59a14f', 'edc948', 'b07aa1', 'ff9da7', '9c755f', 'bab0ac']
    node2colors = {node: [] for node in graph.nodes}
    for j, OG in enumerate(OGs):
        for node in OG:
            node2colors[node].append(cycle[j % 10])

    node_colors = []
    for node in graph.nodes:
        colors = node2colors[node]
        if colors:
            node_colors.append(blend_colors(colors))
        else:
            node_colors.append('#1b1b1b')
    return node_colors


# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = [adj.split(':') for adj in adjs.split(',')]

# Load connected components
CCs = {}
with open('../connect_ggraph/out/gconnect2.txt') as file:
    for line in file:
        CCid, nodes = line.rstrip().split(':')
        CCs[CCid] = set(nodes.split(','))

# Load OGs
OG3s = load_OGs('../subcluster3_ggraph/out/ggraph2/gclusters.txt')
OG4s = load_OGs('../subcluster4_ggraph/out/ggraph2/gclusters.txt')
OG5s = load_OGs('../clique5+_community/out/ggraph2/5clique/gclusters.txt')
OG6s = load_OGs('../clique5+_community/out/ggraph2/6clique/gclusters.txt')

# Make output directory
if not os.path.exists('out/ggraph2/'):
    os.makedirs('out/ggraph2/')  # Recursive folder creation

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
    fig.savefig(f'out/ggraph2/{i}_{CCid}_source.png')
    plt.close()

    # Draw graph labeled by cluster
    for j, OGjs in zip(range(3, 7), [OG3s, OG4s, OG5s, OG6s]):
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
        nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.5)
        nx.draw_networkx_nodes(G.nodes, pos, node_size=node_size, linewidths=0, node_color=get_node_colors(G, OGjs[CCid]))
        fig.tight_layout()
        ax.axis('off')
        fig.savefig(f'out/ggraph2/{i}_{CCid}_OG{j}.png')
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
    fig.savefig(f'out/ggraph2/{i}_{CCid}_edge.png')
    plt.close()

"""
DEPENDENCIES
../clique5+_community/clique5+_community2.py
    ../clique5+_community/out/ggraph2/5clique/gclusters.txt
    ../clique5+_community/out/ggraph2/6clique/gclusters.txt
../connect_ggraph/connect_ggraph2.py
    ../connect_ggraph/out/gconnect2.txt
../hits2ggraph/hits2ggraph2.py
    ../hits2ggraph/out/ggraph2.tsv
../subcluster3_ggraph/subcluster3_ggraph2.py
    ../subcluster3_ggraph/out/ggraph2/gclusters.txt
../subcluster4_ggraph/subcluster4_ggraph2.py
    ../subcluster4_ggraph/out/ggraph2/gclusters.txt
"""