"""Finds connected components of graph, implemented as a depth-first search."""


def get_connected_components(graph):
    # INITIALIZE
    component = None
    components = []
    expand_stack = []  # Stack to expand current component
    search_stack = sorted(graph)  # Stack to search for new component; sort to ensure consistent order
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
