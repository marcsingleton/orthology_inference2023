"""Cluster graph of protein BRHs by triangle criterion, implemented as a variant of a depth-first search."""

from itertools import combinations


class Node:
    def __init__(self, label, adj_labels, marked=False):
        self.label = label
        self.adj_labels = adj_labels
        self.marked = marked

    def __repr__(self):
        return f"Node('{self.label}', {self.adj_labels})"


def OG_member(OG, adj_labels):
    """Return Boolean indicating if adjacency list is member of OG."""

    count = 0
    for adj_label in adj_labels:
        if adj_label in OG:
            count += 1
        if count == 2:
            return True
    return False


def parse_adjacencies(node, OG, graph, expand_stack):
    """Iterate over adjacencies of node, adding node to OG and to expand stack as appropriate."""

    for adj_label in node.adj_labels:
        if adj_label in OG:  # Skip previously added nodes; not strictly necessary but prevents redundant loops
            continue

        adj = graph[adj_label]
        if OG_member(OG, adj.adj_labels):
            OG.add(adj_label)
            if not adj.marked:
                expand_stack.append(adj_label)

    node.marked = True


def cluster(adj_dict):
    # INITIALIZE
    OG = None
    OGs = []
    expand_stack = list()  # Stack to expand current OG
    search_stack = list(adj_dict)  # Stack to search for new OG; initialize with all nodes
    graph = {label: Node(label, adj_labels) for label, adj_labels in adj_dict.items()}  # Convert adj_dict to graph

    # LOOP
    while expand_stack or search_stack:
        # Exhaust expand stack first
        while expand_stack:
            node = graph[expand_stack.pop()]
            parse_adjacencies(node, OG, graph, expand_stack)
        if OG is not None:  # Only record OG if not None; only False in first iteration
            OGs.append(OG)
            OG = None

        # Proceed to search stack when expand stack is empty
        while search_stack and OG is None:
            node = graph[search_stack.pop()]
            if node.marked:  # Skip previously added nodes; necessary to prevent loops and incomplete OGs
                continue

            for adj1_label, adj2_label in combinations(node.adj_labels, 2):  # Check adjacencies of adjacent pairs to create new OG
                adj1, adj2 = graph[adj1_label], graph[adj2_label]
                if adj1.label in adj2.adj_labels:  # Assumes undirected
                    OG = set([node.label, adj1.label, adj2.label])  # Sets check membership in constant time
                    expand_stack.extend([adj1.label, adj2.label])

                    parse_adjacencies(node, OG, graph, expand_stack)

    return OGs
