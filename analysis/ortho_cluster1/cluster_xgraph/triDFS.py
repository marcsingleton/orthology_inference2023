"""Cluster graph of protein BRHs by triangle criterion, implemented as a variant of a depth-first search."""

from itertools import combinations


def cluster(graph):
    # INITIALIZE
    OG = set()
    OGs = []
    expand_stack = []  # Stack to expand current OG
    search_stack = sorted([(node, adj) for node, adjs in graph.items() for adj in adjs])  # Stack to search for new OG; sort to ensure consistent order
    marked = set()

    # LOOP
    while expand_stack or search_stack:
        # Exhaust expand stack first
        while expand_stack:
            edge = expand_stack.pop()
            if edge in marked:  # Prevents infinite loops
                continue

            node1, node2 = edge
            for node3 in graph[node1]:
                if node3 in graph[node2]:  # Assumes undirected
                    edges = [frozenset([node1, node3]), frozenset([node2, node3])]
                    OG |= set(edges)  # Sets check membership in constant time
                    expand_stack.extend(edges)
            marked.add(edge)
        if OG:  # Only record OG has members; only False in first iteration to prevent adding "empty" OG
            OGs.append(OG)
            OG = set()

        # Proceed to search stack when expand stack is empty
        while search_stack and not OG:
            edge = search_stack.pop()
            if frozenset(edge) in marked:  # Prevents revisiting previous OGs
                continue

            node1, node2 = edge
            for node3 in sorted(graph[node1]):  # Sort to ensure consistent order
                if node3 in graph[node2]:  # Assumes undirected
                    edges = [frozenset([label1, label2]) for label1, label2 in combinations([node1, node2, node3], 2)]
                    OG |= set(edges)  # Sets check membership in constant time
                    expand_stack.extend(edges)
                    break

    return OGs
