"""Cluster graph of gene BRHs by triangle criterion, implemented as a variant of a depth-first search."""

from itertools import combinations


class Node:
    def __init__(self, label, adj_gns, marked=False):
        self.label = label
        self.adj_gns = adj_gns
        self.marked = marked

    def __repr__(self):
        return f"Node('{self.label}', {list(self.adj_gns.keys())})"


def get_OG_links(OG, adj_gns):
    """Return list of edges to current OG from an adjacency dictionary."""

    links = []
    for adj_label in adj_gns:
        if adj_label in OG:
            links.append(adj_label)

    return links


def parse_adjacencies(node, OG, graph, expand_stack):
    """Iterate over adjacencies of node, adding them to OG and expand stack as appropriate."""

    for adj_label in node.adj_gns:
        if adj_label in OG:  # Skip previously added nodes; not strictly necessary but prevents redundant loops
            continue

        adj = graph[adj_label]
        links = get_OG_links(OG, adj.adj_gns)
        if len(links) >= 2:
            OG[adj.label] = {}
            for link in links:
                OG[adj.label][link] = adj.adj_gns[link]
                OG[link][adj.label] = graph[link].adj_gns[adj.label]
            if not adj.marked:
                expand_stack.append(adj_label)

    node.marked = True


def cluster(adj_dict):
    # INITIALIZE
    OG = None
    OGs = []
    expand_stack = list()  # Stack to expand current OG
    search_stack = list(adj_dict)  # Stack to search for new OG; initialize with all nodes
    graph = {label: Node(label, adj_gns) for label, adj_gns in adj_dict.items()}  # Convert adj_dict to graph

    # LOOP
    while expand_stack or search_stack:
        # Exhaust expand stack first
        while expand_stack:
            node = graph[expand_stack.pop()]
            parse_adjacencies(node, OG, graph, expand_stack)
        if OG is not None:  # Only record OG if not None; only False in first iteration to prevent adding "empty" OG
            OGs.append(OG)
            OG = None

        # Proceed to search stack when expand stack is empty
        while search_stack and OG is None:
            node = graph[search_stack.pop()]
            if node.marked:  # Skip previously added nodes; necessary to prevent loops and incomplete OGs
                continue

            for adj1_label, adj2_label in combinations(node.adj_gns, 2):  # Check adjacencies of adjacent pairs to create new OG
                adj1, adj2 = graph[adj1_label], graph[adj2_label]
                if adj1.label in adj2.adj_gns:  # Assumes undirected
                    OG = {node.label: {adj1.label: node.adj_gns[adj1.label], adj2.label: node.adj_gns[adj2.label]},
                          adj1.label: {node.label: adj1.adj_gns[node.label], adj2.label: adj1.adj_gns[adj2.label]},
                          adj2.label: {adj1.label: adj2.adj_gns[adj1.label], node.label: adj2.adj_gns[node.label]}}
                    expand_stack.extend([adj1.label, adj2.label])

                    parse_adjacencies(node, OG, graph, expand_stack)

    return OGs
