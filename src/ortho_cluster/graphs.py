"""Functions for finding structures in graphs."""

from itertools import combinations


def get_connected_components(graph):
    """Find connected components of graph, implemented as a depth-first search."""
    # INITIALIZE
    component = None
    components = []
    expand_stack = []  # Stack to expand current component
    search_stack = list(graph)  # Stack to search for new component
    marked = set()

    # LOOP
    while expand_stack or search_stack or component is not None:
        # Exhaust expand stack first
        while expand_stack:
            node = expand_stack.pop()
            if node in marked:
                continue
            component.add(node)
            expand_stack.extend(graph[node])
            marked.add(node)
        if component is not None:  # Only record component if not None; only False in first iteration
            components.append(component)
            component = None

        # Proceed to search stack when expand stack is empty
        while search_stack and component is None:
            node = search_stack.pop()
            if node in marked:  # Skip previously added nodes; necessary to prevent loops and incomplete components
                continue
            component = {node}
            expand_stack.extend(graph[node])
    return components


def get_triangle_clusters(graph):
    """Cluster graph of protein BRHs by triangle criterion, implemented as a variant of a depth-first search."""
    # INITIALIZE
    cluster = set()
    clusters = []
    expand_stack = []  # Stack to expand current cluster
    search_stack = [(node, adj) for node, adjs in graph.items() for adj in adjs]  # Stack to search for new cluster;
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
                    cluster |= set(edges)  # Sets check membership in constant time
                    expand_stack.extend(edges)
            marked.add(edge)
        if cluster:  # Only record cluster has members; only False in first iteration to prevent adding "empty" cluster
            clusters.append(cluster)
            cluster = set()

        # Proceed to search stack when expand stack is empty
        while search_stack and not cluster:
            edge = search_stack.pop()
            if frozenset(edge) in marked:  # Prevents revisiting previous OGs
                continue

            node1, node2 = edge
            for node3 in sorted(graph[node1]):  # Sort to ensure consistent order
                if node3 in graph[node2]:  # Assumes undirected
                    edges = [frozenset([label1, label2]) for label1, label2 in combinations([node1, node2, node3], 2)]
                    cluster |= set(edges)  # Sets check membership in constant time
                    expand_stack.extend(edges)
                    break

    return clusters
