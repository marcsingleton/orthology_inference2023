"""Cluster ggraphs by 4-clique percolation."""

import matplotlib.pyplot as plt
import os
from itertools import combinations


def get_triangles(graph):
    triangles = set()
    for node, adjs in graph.items():
        for adj1, adj2 in combinations(adjs, 2):
            if adj1 in graph[adj2] and adj2 in graph[adj1]:
                triangles.add(frozenset([node, adj1, adj2]))
    return triangles


def get_edges(cluster):
    edges = set()
    for triangle in cluster:
        edges.update(set(combinations(triangle, 2)))
    return edges


def cluster(graph):
    # INITIALIZE
    OG = set()
    OGs = []
    expand_stack = set()  # Stack to expand current OG
    search_stack = sorted(get_triangles(graph))  # Stack to search for new OG; sort to ensure consistent order
    marked = set()

    # LOOP
    while expand_stack or search_stack:
        # Exhaust expand stack first
        while expand_stack:
            triangle = expand_stack.pop()
            if triangle in marked:  # Prevents infinite loops
                continue

            node1, node2, node3 = triangle
            for node4 in (graph[node1] & graph[node2] & graph[node3]):
                triangles = set([frozenset([node4, nodeA, nodeB]) for nodeA, nodeB in combinations(triangle, 2)])
                OG |= triangles
                expand_stack |= triangles
            marked.add(triangle)

        if OG:  # Only record OG has members; only False in first iteration to prevent adding "empty" OG
            OGs.append(OG)
            OG = set()

        # Proceed to search stack when expand stack is empty
        while search_stack and not OG:
            triangle = search_stack.pop()
            if frozenset(triangle) in marked:  # Prevents revisiting previous OGs
                continue

            node1, node2, node3 = triangle
            for node4 in sorted(graph[node1]):  # Sort to ensure consistent order
                if node4 in graph[node2] and node4 in graph[node3]:  # Assumes undirected
                    triangles = set([triangle, *[frozenset([node4, nodeA, nodeB]) for nodeA, nodeB in combinations(triangle, 2)]])
                    OG |= triangles  # Sets check membership in constant time
                    expand_stack |= triangles
                    break

    return OGs


# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = set([adj.split(':')[0] for adj in adjs.split(',')])

# Load connected components
CCs = []
with open('../connect_ggraph/out/gconnect2.txt') as file:
    for line in file:
        _, nodes = line.rstrip().split(':')
        CCs.append(set(nodes.split(',')))

OGs = []
CCtypes = [{} for _ in range(5)]
for CC in CCs:
    subggraph = {node: ggraph[node] for node in CC}

    # Cluster by triangle criterion
    subOGs = [get_edges(c) for c in cluster(subggraph)]
    OGs.append(subOGs)

    # Classify CCs
    subnOGs = [set([node for edge in subOG for node in edge]) for subOG in subOGs]
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

# Make plots output directory
if not os.path.exists('out/ggraph2/'):
    os.makedirs('out/ggraph2/')  # Recursive folder creation

# Write OGs to file
j = 0
with open('out/ggraph2/gclusters.txt', 'w') as file:
    for i, subOGs in enumerate(OGs):
        CCid = hex(i)[2:].zfill(4)
        for OG in subOGs:
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
plt.savefig('out/ggraph2/connectnum-OGnum_type_dist1-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/ggraph2/connectnum-OGnum_type_dist1-2.png')
plt.close()

plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3', color='C3')
plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4', color='C4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components across number of OGs')
plt.legend()
plt.savefig('out/ggraph2/connectnum-OGnum_type_dist2-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/ggraph2/connectnum-OGnum_type_dist2-2.png')
plt.close()

plt.pie([sum(CCtype.values()) for CCtype in CCtypes], labels=[f'Type {i}' for i in range(len(CCtypes))])
plt.title('Connected components by type')
plt.savefig('out/ggraph2/type_pie.png')
plt.close()

for i, CCtype in enumerate(CCtypes):
    print(f'Type {i}:', sum(CCtype.values()))

"""
OUTPUT
Type 0: 1319
Type 1: 11080
Type 2: 1838
Type 3: 232
Type 4: 691

DEPENDENCIES
../connect_ggraph/connect_ggraph.py
    ../connect_ggraph/out/gconnect2.txt
../hits2ggraph/hits2ggraph.py
    ../hits2ggraph/out/ggraph2.tsv
"""