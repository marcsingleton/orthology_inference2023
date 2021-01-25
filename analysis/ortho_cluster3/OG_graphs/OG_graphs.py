"""Draw largest OGs."""

import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import os
from math import exp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace


def load_pOGs(path):
    pOGs = {}
    with open(path) as file:
        for line in file:
            OGid, _, edges = line.rstrip().split(':')
            gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
            try:
                pOGs[OGid].append(gnids)
            except KeyError:
                pOGs[OGid] = [gnids]
    return pOGs


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


# Load seq metadata
gnid2ppids = {}
ppid2gnid = {}
ppid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, repr = line.split()
        ppid2spid[ppid] = spid
        ppid2gnid[ppid] = gnid
        if repr == 'True':
            try:
                gnid2ppids[gnid].append(ppid)
            except KeyError:
                gnid2ppids[gnid] = [ppid]

# Load pgraph
pgraph = {}
with open('../hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        pgraph[node] = bitscores

# Load OGs
OGs = {}
with open('../clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

# Load pOGs
pOG3s = load_pOGs('../subcluster_pgraph/out/pclusters.txt')
pOG4s = load_pOGs('../clique4+_pcommunity/out/4clique/pclusters.txt')
pOG5s = load_pOGs('../clique4+_pcommunity/out/5clique/pclusters.txt')
pOG6s = load_pOGs('../clique4+_pcommunity/out/6clique/pclusters.txt')

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = sorted(OGs, key=lambda x: len(OGs[x]), reverse=True)
for i, OGid in enumerate(OGids[:50]):  # 50 largest OGs
    # Make subpgraph
    subpgraph = {}
    ppids = set([ppid for gnid in OGs[OGid] for ppid in gnid2ppids[gnid]])
    for ppid in ppids:
        # Collect hits by SPID
        spids = {}
        for adj_node, adj_bitscore in pgraph.get(ppid, {}).items():  # In case PPID has no hits
            spid = ppid2spid[adj_node]
            try:
                spids[spid][adj_node] = adj_bitscore
            except KeyError:
                spids[spid] = {adj_node: adj_bitscore}

        # Find max hits for each SPID
        d = dict()
        for spid, adjs in spids.items():
            max_bitscore = max(adjs.values())
            d.update({adj_node: adj_bitscore for adj_node, adj_bitscore in adjs.items() if adj_bitscore == max_bitscore})
        subpgraph[ppid] = d

    # Filter subpgraph by reciprocity
    for node, adjs in subpgraph.items():
        del_keys = []
        for adj in adjs:
            if not (adj in subpgraph and node in subpgraph[adj]):
                del_keys.append(adj)
        for del_key in del_keys:
            del adjs[del_key]

    # Create graph
    G = nx.Graph()
    for node, adjs in subpgraph.items():
        G.add_node(node)
        for adj, w in adjs.items():
            if (node, adj) in G.edges:
                edge_data = G.get_edge_data(node, adj)
                edge_data['weight'] = edge_data.get('weight', 0) + float(w)  # Sum edge weights
            else:
                G.add_edge(node, adj, weight=float(w))

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
    node_size = 25/(1 + exp(0.01*(len(subpgraph)-400))) + 10  # Adjust node size
    FB = [node for node in G.nodes if node.startswith('FBpp')]
    NCBI = G.nodes - FB

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.5)
    nx.draw_networkx_nodes(NCBI, pos, node_size=node_size, linewidths=0, node_color='C0', label='NCBI')
    nx.draw_networkx_nodes(FB, pos, node_size=node_size, linewidths=0, node_color='C1', label='FlyBase')

    fig.legend(markerscale=(1 if node_size > 22.5 else 22.5/node_size), loc=locs[-1])
    fig.tight_layout()
    ax.axis('off')
    fig.savefig(f'out/{i}_{OGid}_source.png')
    plt.close()

    # Draw graph labeled by cluster
    for j, OGjs in zip(range(3, 7), [pOG3s, pOG4s, pOG5s, pOG6s]):
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
        nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.5)
        nx.draw_networkx_nodes(G.nodes, pos, node_size=node_size, linewidths=0, node_color=get_node_colors(G, OGjs[OGid]))
        fig.tight_layout()
        ax.axis('off')
        fig.savefig(f'out/{i}_{OGid}_OG{j}.png')
        plt.close()

    # Draw graph labeled by edge
    node_size = 20/(1 + exp(0.01*(len(subpgraph)-400))) + 10  # Adjust node size
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
    fig.savefig(f'out/{i}_{OGid}_edge.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_search/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_pcommunity/clique4+_pcommunity.py
    ../clique4+_pcommunity/out/4clique/pclusters.txt
    ../clique4+_pcommunity/out/5clique/pclusters.txt
    ../clique4+_pcommunity/out/6clique/pclusters.txt
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph2.tsv
"""