"""Clusters graph by k-clique percolation."""

from collections import defaultdict
from itertools import combinations

import networkx as nx


# Time out errors and handlers
class CliqueError(Exception):
    pass


class PercolateError(Exception):
    pass


def clique_handler(signum, frame):
    raise CliqueError


def percolate_handler(signum, frame):
    raise PercolateError


# Graph functions
def k_clique_communities(G, k, cliques=None):
    """Find k-clique communities in graph using the percolation method.

    A k-clique community is the union of all cliques of size k that
    can be reached through adjacent (sharing k-1 nodes) k-cliques.

    This function is taken nearly verbatim from networkx with the minor
    change that it reports edges rather than nodes.

    Parameters
    ----------
    G : NetworkX graph

    k : int
       Size of smallest clique

    cliques: list or generator
       Precomputed cliques (use networkx.find_cliques(G))

    Returns
    -------
    Yields sets of nodes, one for each k-clique community.

    References
    ----------
    .. [1] Gergely Palla, Imre Derényi, Illés Farkas1, and Tamás Vicsek,
       Uncovering the overlapping community structure of complex networks
       in nature and society Nature 435, 814-818, 2005,
       doi:10.1038/nature03607
    """
    if k < 2:
        raise nx.NetworkXError(f'k={k}, k must be greater than 1.')
    if cliques is None:
        cliques = nx.find_cliques(G)
    cliques = [frozenset(c) for c in cliques if len(c) >= k]

    # First index which nodes are in which cliques
    membership_dict = defaultdict(list)
    for clique in cliques:
        for node in clique:
            membership_dict[node].append(clique)

    # For each clique, see which adjacent cliques percolate
    perc_graph = nx.Graph()
    perc_graph.add_nodes_from(cliques)
    for clique in cliques:
        for adj_clique in _get_adjacent_cliques(clique, membership_dict):
            if len(clique.intersection(adj_clique)) >= (k - 1):
                perc_graph.add_edge(clique, adj_clique)

    # Connected components of clique graph with perc edges
    # are the percolated cliques
    for component in nx.connected_components(perc_graph):
        yield frozenset([frozenset(edge) for clique in component for edge in combinations(clique, 2)])


def _get_adjacent_cliques(clique, membership_dict):
    adjacent_cliques = set()
    for n in clique:
        for adj_clique in membership_dict[n]:
            if clique != adj_clique:
                adjacent_cliques.add(adj_clique)
    return adjacent_cliques


def k_clique_communities_progressive(G, k, cliques=None):
    """Find k-clique communities in graph using the percolation method.

    This method differs from the networkx version in two respects. The first is
    it returns edges rather than nodes. The second, and more important
    difference is that it does not exhaustively find edges in the percolation
    graph. Since in some cases all cliques are found relatively quickly, the
    bottleneck is calculating the C*(C-1)/2 overlaps between all cliques,
    particularly if the number of cliques, C, is large. Since joining a
    k-clique community only requires that a clique is connected to a community,
    this approach is needlessly expensive for large graphs.

    This method assumes the number of connected components is small and
    therefore progressively checks each subsequent clique against a list of
    currently known connected components. All components which have at least
    one clique with sufficient overlap with the current clique are merged. The
    search can start with an arbitrary clique as the first connected component.
    This implementation, however, first sorts them by the number of edges
    associated with each node in the clique which are also not a member of that
    clique. This begins the search with cliques whose nodes are highly
    connected to nodes outside of that clique and are thus likely to have
    sufficient overlap with other cliques. The order of cliques within
    components reflects the history of merges and thus will contain runs of
    sorted cliques but will not be sorted globally. Using a priority queue to
    ensure highly-connected cliques are checked first within each component
    could further increase performance.

    Additionally, cliques without overlap are inserted at the end of the list
    of connected components. This ensures highly connected components remain at
    the front of the list and are thus checked first.

    At worst, the algorithm checks every clique against every other clique,
    which is the same time complexity as the naive approach. (It does not store
    the edges of the percolation graph, however, yielding an advantage in terms
    of space complexity.) At best, there is one connected component and the
    cliques are sorted so the initial seed clique overlaps with all others. In
    this case the complexity is linear in the number of cliques.
    """
    if k < 2:
        raise nx.NetworkXError(f'k={k}, k must be greater than 1.')
    if cliques is None:
        cliques = nx.find_cliques(G)
    cliques = [frozenset(c) for c in cliques if len(c) >= k]
    cliques = sorted(cliques, key=lambda c: sum([len(G[node]) for node in c]) - len(c)*(len(c)-1))

    components = [[cliques.pop()]] if cliques else []  # List of lists of cliques (frozensets of nodes)
    for clique in cliques:
        intersect = [clique]
        disjoint = []
        for component in components:
            if k_percolates(clique, component, k):
                intersect.extend(component)
            else:
                disjoint.append(component)
        components = [intersect]
        components.extend(disjoint)

    for component in components:
        yield frozenset([frozenset(edge) for clique in component for edge in combinations(clique, 2)])


def k_percolates(clique, CC, k):
    for component in CC:
        if len(clique & component) >= k - 1:
            return True
    return False
