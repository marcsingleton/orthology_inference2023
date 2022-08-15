"""Draw largest connected components."""

import os
from math import floor, log10

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from numpy import linspace


class ColorCycler:
    def __init__(self, color_cycle):
        self.color_cycle = color_cycle.copy()
        self.color_cycle_template = color_cycle.copy()

    def get_color(self, exclude_color=None):
        if len(self.color_cycle) == 0:
            self.color_cycle = self.color_cycle_template.copy()
        elif len(self.color_cycle) == 1 and exclude_color == self.color_cycle[0]:
            self.color_cycle = self.color_cycle_template.copy()

        idx = 0
        for color in self.color_cycle:
            if color != exclude_color:
                break
            idx += 1

        return self.color_cycle.pop(idx)

    def add_color(self, color):
        self.color_cycle.insert(0, color)


def load_OGs(path):
    OGs = {}
    with open(path) as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            component_id, OGid = fields['component_id'], fields['OGid']
            ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
            try:
                OGs[component_id][OGid] = ppids
            except KeyError:
                OGs[component_id] = {OGid: ppids}
    return OGs


def blend_colors(colors):
    colors = [color.lstrip('#') for color in colors]
    rgbs = [[int(color[i:i+2], 16) for i in range(0, 6, 2)] for color in colors]
    avg = [int(sum(c) / len(c)) for c in zip(*rgbs)]
    return '#' + ''.join([f'{c:X}' for c in avg])


def get_node_colors(OGs_k, graph, color_cycler):
    # 1 Assign colors to OGs
    OG_colors_k = []

    # 1.1 Get colors for smallest k by color cycle
    OGs0 = OGs_k[0]
    OG_colors0 = {}
    for OGid, OG in sorted(OGs0.items(), key=lambda x: -len(x[1])):
        color = color_cycler.get_color()
        OG_colors0[OGid] = color
    OG_colors_k.append(OG_colors0)

    for OGs in OGs_k[1:]:
        # 1.2 Find parent OGs
        OGid02OGids = {}
        for OGid, OG in OGs.items():
            for OGid0, OG0 in OGs0.items():
                if set(OG0) >= set(OG):
                    try:
                        OGid02OGids[OGid0].append(OGid)
                    except KeyError:
                        OGid02OGids[OGid0] = [OGid]
                    break

        # 1.3 Return colors of OGs without children to color cycle
        for OGid0, color in OG_colors0.items():
            if OGid0 not in OGid02OGids:
                color_cycler.add_color(color)

        # 1.4 Assign colors
        OG_colors = {}
        for OGid0, OGids in OGid02OGids.items():
            heir_OGid = max(OGids, key=lambda x: len(OGs[x]))
            heir_color = OG_colors0[OGid0]
            for OGid in OGids:
                if OGid == heir_OGid:
                    OG_colors[OGid] = heir_color  # Largest OG receives parent color
                else:
                    color = color_cycler.get_color(exclude_color=heir_color)  # New OGs receive new color
                    OG_colors[OGid] = color
        for OGid in OGs:  # Catch any remaining OGs (should only occur in cases where clique percolation timed out)
            if OGid not in OG_colors:
                color = color_cycler.get_color()  # New OGs receive new color
                OG_colors[OGid] = color
        OG_colors_k.append(OG_colors)
        OGs0 = OGs
        OG_colors0 = OG_colors

    # 2 Color nodes
    node_colors_k = []
    for OGs, OG_colors in zip(OGs_k, OG_colors_k):
        # 2.1 Collect node colors
        node2colors = {node: [] for node in graph.nodes}
        for OGid, color in OG_colors.items():
            for node in OGs[OGid]:
                node2colors[node].append(color)

        # 2.2 Blend node colors
        node_colors = []
        for node in graph.nodes:
            colors = node2colors[node]
            if colors:
                node_colors.append(blend_colors(colors))
            else:
                node_colors.append(null_color)
        node_colors_k.append(node_colors)

    return node_colors_k


num_components = 100

fig_width = 7.5
dpi = 400

margin = 0
margin_data = 0.1
axes_rect = (margin, margin, 1 - margin, 1 - margin)

color_cycle = ['#4E79A7F8', '#F28E2BF8', '#E15759F8', '#499894F8', '#59A14FF8', '#B6992DF8', '#B07AA1F8', '#D37295F8',
               '#A0CBE8F8', '#FFBE7DF8', '#FF9D9AF8', '#86BCB6F8', '#8CD17DF8', '#F1CE63F8', '#D4A6C8F8', '#FABFD2F8']
null_color = '#BFBFBF'
node_slope1 = 175
node_intercept1 = 0
edge_alpha1 = 0.5
edge_width1 = 0.25

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

# Load OGs
OGs_3 = load_OGs('../cluster3_graph/out/clusters.tsv')
OGs_4 = load_OGs('../cluster4+_graph/out/4clique/clusters.tsv')
OGs_5 = load_OGs('../cluster4+_graph/out/5clique/clusters.tsv')
OGs_6 = load_OGs('../cluster4+_graph/out/6clique/clusters.tsv')

ks = list(range(3, 7))
OGs_k = [OGs_3, OGs_4, OGs_5, OGs_6]

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
    node_size1 = node_slope1 * fig_width * fig_height / len(subgraph) + node_intercept1  # Node size is inversely proportional to node density
    node_size2 = node_slope2 * fig_width * fig_height / len(subgraph) + node_intercept2

    # Draw graph labeled by cluster
    node_colors_k = get_node_colors([OGs[component_id] for OGs in OGs_k], nx_graph, ColorCycler(color_cycle))
    for k, node_colors in zip(ks, node_colors_k):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes(axes_rect)

        nx.draw_networkx_edges(nx_graph, positions, ax=ax, edge_color=null_color, alpha=edge_alpha1, width=edge_width1)
        nx.draw_networkx_nodes(nx_graph, positions, ax=ax, node_size=node_size1, linewidths=0, node_color=node_colors)
        ax.set_xlim((xmin - margin_data * xlen, xmax + margin_data * xlen))  # Set manually because draw_networkx_edges hard codes the data limits with 5% padding
        ax.set_ylim((ymin - margin_data * ylen, ymax + margin_data * ylen))

        ax.axis('off')
        fig.savefig(f'out/{i:02}_{component_id}_OG{k}.png', dpi=dpi)
        plt.close()

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
../cluster3_graph/cluster.py
    ../cluster3_graph/out/clusters.tsv
../cluster4+_graph/cluster.py
    ../cluster4+_graph/out/4clique/clusters.tsv
    ../cluster4+_graph/out/5clique/clusters.tsv
    ../cluster4+_graph/out/6clique/clusters.tsv
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""