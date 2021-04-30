"""Finds connected components of graph, implemented as a depth-first search."""


class Node:
    def __init__(self, label, adj_labels, marked=False):
        self.label = label
        self.adj_labels = adj_labels
        self.marked = marked

    def __repr__(self):
        return f"Node('{self.label}', {self.adj_labels})"


def connect(adj_dict):
    # INITIALIZE
    component = None
    components = []
    expand_stack = list()  # Stack to expand current component
    search_stack = sorted(list(adj_dict))  # Stack to search for new component; sort to ensure consistent order
    graph = {label: Node(label, adj_labels) for label, adj_labels in adj_dict.items()}  # Convert adj_dict to graph

    # LOOP
    while expand_stack or search_stack or component is not None:
        # Exhaust expand stack first
        while expand_stack:
            node = graph[expand_stack.pop()]
            if node.marked:
                continue
            component.add(node.label)
            expand_stack.extend(node.adj_labels)
            node.marked = True
        if component is not None:  # Only record component if not None; only False in first iteration
            components.append(component)
            component = None

        # Proceed to search stack when expand stack is empty
        while search_stack and component is None:
            node = graph[search_stack.pop()]
            if node.marked:  # Skip previously added nodes; necessary to prevent loops and incomplete components
                continue
            component = set([node.label])
            expand_stack.extend(node.adj_labels)
    return components
