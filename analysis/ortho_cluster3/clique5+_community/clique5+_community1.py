"""Cluster by k-clique percolation."""

import os
import signal
from collections import defaultdict
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.components import connected_components
from networkx.algorithms.core import k_core


# Define time out errors and handlers
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
    it returns reports rather than nodes. The second, and more important
    difference is that it does not exhaustively find edges in the percolation
    graph. Since in some cases all the cliques are found relatively quickly,
    the bottleneck is calculating the C*(C-1)/2 overlaps between all cliques,
    particularly if the number of cliques, C, is large. Since joining a
    k-clique community only requires that a clique is connected to a community,
    this approach is needlessly expensive for large graphs.

    This method assumes the number of connected components is small and
    therefore progressively constructs by checking each subsequent clique
    against a list of currently known connected components. All components
    which have at least one clique with sufficient overlap with the current
    clique are merged. The search can start with an arbitrary clique as the
    first connected component. This implementation, however, first sorts them
    by the number of edges associated with each node in the clique which are
    also not a member of that clique. This begins the search with cliques whose
    nodes are highly connected to nodes outside of that clique and are thus
    likely to have sufficient overlap with other cliques. The order of cliques
    within components reflects the history of merges and thus will contain runs
    of sorted cliques but will not be sorted globally. Using a priority queue
    to ensure highly-connected cliques are checked first within each component
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

    CCs = [[cliques.pop()]]  # List of lists of cliques (frozensets of nodes)
    for clique in cliques:
        intersect = [clique]
        disjoint = []
        for CC in CCs:
            if k_percolates(clique, CC, k):
                intersect.extend(CC)
            else:
                disjoint.append(CC)
        CCs = [intersect]
        CCs.extend(disjoint)

    for component in CCs:
        yield frozenset([frozenset(edge) for clique in component for edge in combinations(clique, 2)])


def edges2nodes(edges):
    return set([node for edge in edges for node in edge])


def k_percolates(clique, CC, k):
    for component in CC:
        if len(clique & component) > k - 1:
            return True
    return False


# Output functions
def classify_CC(CCtypes, subOGs):
    subnOGs = [edges2nodes(subOG) for subOG in subOGs]
    if len(subnOGs) == 0:
        CCtypes[0][len(subnOGs)] = CCtypes[0].get(len(subnOGs), 0) + 1  # Component has no OGs
    elif len(subnOGs) == 1:
        if len(subnOGs[0]) == len(CC):
            CCtypes[1][len(subnOGs)] = CCtypes[1].get(len(subnOGs), 0) + 1  # Component and OG are equal
        else:
            CCtypes[2][len(subnOGs)] = CCtypes[2].get(len(subnOGs), 0) + 1  # Component has single OG which is a subset of the component
    elif any([set.intersection(nOG1, nOG2) for nOG1, nOG2 in combinations(subnOGs, 2)]):
        CCtypes[4][len(subnOGs)] = CCtypes[4].get(len(subnOGs), 0) + 1  # Component has multiple non-disjoint OGs
    else:
        CCtypes[3][len(subnOGs)] = CCtypes[3].get(len(subnOGs), 0) + 1  # Component has multiple pairwise disjoint OGs


def save_results(OGs, CCtypes, k):
    # Make plots output directory
    if not os.path.exists(f'out/ggraph1/{k}clique/'):
        os.makedirs(f'out/ggraph1/{k}clique/')  # Recursive folder creation

    # Write OGs to file
    j = 0
    with open(f'out/ggraph1/{k}clique/gclusters.txt', 'w') as file:
        for i, subOGs in enumerate(OGs):
            CCid = hex(i)[2:].zfill(4)
            for OG in sorted(subOGs, key=lambda x: sorted(edges2nodes(x))):  # Ensure consistent ids by keying on sorted node list
                OGid = hex(j)[2:].zfill(4)
                file.write(CCid + ':' + OGid + ':' + '\t'.join([f'{node1},{node2}' for node1, node2 in OG]) + '\n')
                j += 1

    # Plots
    plt.bar(CCtypes[0].keys(), CCtypes[0].values(), label='Type 0')
    plt.bar(CCtypes[1].keys(), CCtypes[1].values(), label='Type 1')
    plt.bar(CCtypes[2].keys(), CCtypes[2].values(), bottom=CCtypes[1][1], label='Type 2')
    plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3')
    plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4')
    plt.xlabel('Number of OGs in connected component')
    plt.ylabel('Number of connected components')
    plt.title('Distribution of connected components across number of OGs')
    plt.legend()
    plt.savefig(F'out/ggraph1/{k}clique/connectnum-OGnum_type_dist1-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/ggraph1/{k}clique/connectnum-OGnum_type_dist1-2.png')
    plt.close()

    plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3', color='C3')
    plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4', color='C4')
    plt.xlabel('Number of OGs in connected component')
    plt.ylabel('Number of connected components')
    plt.title('Distribution of connected components across number of OGs')
    plt.legend()
    plt.savefig(f'out/ggraph1/{k}clique//connectnum-OGnum_type_dist2-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/ggraph1/{k}clique/connectnum-OGnum_type_dist2-2.png')
    plt.close()

    plt.pie([sum(CCtype.values()) for CCtype in CCtypes], labels=[f'Type {i}' for i in range(len(CCtypes))])
    plt.title('Connected components by type')
    plt.savefig(f'out/ggraph1/{k}clique/type_pie.png')
    plt.close()

    print()
    print(f'{k}-CLIQUE')
    for i, CCtype in enumerate(CCtypes):
        print(f'Type {i}:', sum(CCtype.values()))


# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph1.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = set([adj.split(':')[0] for adj in adjs.split(',')])

# Load connected components
CCs = []
with open('../connect_ggraph/out/gconnect1.txt') as file:
    for line in file:
        CCid, nodes = line.rstrip().split(':')
        CCs.append((CCid, set(nodes.split(','))))

# Load OG4s
OG4s = {}
with open('../subcluster4_ggraph/out/ggraph1/gclusters.txt') as file:
    for line in file:
        CCid, _, edges = line.rstrip().split(':')
        try:
            OG4s[CCid].append([edge.split(',') for edge in edges.split('\t')])
        except KeyError:
            OG4s[CCid] = [[edge.split(',') for edge in edges.split('\t')]]

OG5s = []
OG6s = []
CCtypes5 = [{} for _ in range(5)]
CCtypes6 = [{} for _ in range(5)]
for CCid, CC in CCs:
    # Create graph
    G = nx.Graph()
    for node in CC:
        G.add_node(node)
        for adj in ggraph[node]:
            G.add_edge(node, adj)

    # Cluster by clique percolation
    try:
        # Handle cliques
        signal.signal(signal.SIGALRM, clique_handler)
        if CCid in ['0869', '08a6', '08a7', '08a8', '08a9']:
            raise CliqueError
        signal.alarm(30)
        cliques = list(find_cliques(G))
        signal.alarm(0)

        # Handle percolation
        signal.signal(signal.SIGALRM, percolate_handler)
        signal.alarm(30)
        subOG5s = list(k_clique_communities(G, 5, cliques))
        subOG6s = list(k_clique_communities(G, 6, cliques))
        signal.alarm(0)
    except CliqueError:
        print('CliqueError:', CCid)
        subOG5s = set()
        subOG6s = set()
        for OG4 in OG4s[CCid]:
            G = nx.Graph(OG4)

            core5 = k_core(G, 5)
            for component in connected_components(core5):
                subOG5s.add(frozenset([frozenset(edge) for edge in core5.edges(component)]))
            core6 = k_core(G, 6)
            for component in connected_components(core6):
                subOG6s.add(frozenset([frozenset(edge) for edge in core6.edges(component)]))
    except PercolateError:
        print('PercolateError:', CCid)
        subOG5s = list(k_clique_communities_progressive(G, 5, cliques))
        subOG6s = list(k_clique_communities_progressive(G, 6, cliques))
    OG5s.append(subOG5s)
    OG6s.append(subOG6s)

    # Classify CCs
    classify_CC(CCtypes5, subOG5s)
    classify_CC(CCtypes6, subOG6s)

save_results(OG5s, CCtypes5, 5)
save_results(OG6s, CCtypes6, 6)

"""
OUTPUT
CliqueError: 0869
CliqueError: 08a6
CliqueError: 08a7
CliqueError: 08a8
CliqueError: 08a9
PercolateError: 0cf6
PercolateError: 2ff2

5-CLIQUE
Type 0: 1677
Type 1: 10509
Type 2: 1444
Type 3: 243
Type 4: 591

6-CLIQUE
Type 0: 2000
Type 1: 10190
Type 2: 1478
Type 3: 247
Type 4: 549

DEPENDENCIES
../connect_ggraph/connect_ggraph.py
    ../connect_ggraph/out/gconnect1.txt
../hits2ggraph/hits2ggraph.py
    ../hits2ggraph/out/ggraph1.tsv
"""