"""Plot statistics related to OGs and connected components."""

from itertools import combinations

import matplotlib.pyplot as plt
import pandas as pd

# Load connected components
components = []
with open('../connect_hit_graph/out/components.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        components.append((fields['component_id'], set(fields['ppids'].split(','))))

# Load OGs
component2OGs = {}
with open('out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        component_id = fields['component_id']
        edges = [edge.split(':') for edge in fields['edges'].split(',')]
        try:
            component2OGs[component_id].append(edges)
        except KeyError:
            component2OGs[component_id] = [edges]

# Classify components
rows = []
for component_id, component in components:
    OGs = component2OGs.get(component_id, [])  # Not all components have OGs
    node_sets = [{node for edge in OG for node in edge} for OG in OGs]
    if len(node_sets) == 0:
        component_type = 0  # Component has no OGs
    elif len(node_sets) == 1:
        if len(node_sets[0]) == len(component):
            component_type = 1  # Component and OG are equal
        else:
            component_type = 2  # Component has single OG which is a subset of the component
    elif any([set.intersection(set1, set2) for set1, set2 in combinations(node_sets, 2)]):
        component_type = 4  # Component has multiple non-disjoint OGs
    else:
        component_type = 3  # Component has multiple pairwise disjoint OGs

    rows.append({'component_id': component_id, 'component_type': component_type, 'OGnum': len(OGs)})

# Plots
df = pd.DataFrame(rows)
component_types = [df.loc[df['component_type'] == i, 'OGnum'].value_counts() for i in range(5)]

plt.bar(component_types[0].index, component_types[0].values, label='Type 0')
plt.bar(component_types[1].index, component_types[1].values, label='Type 1')
plt.bar(component_types[2].index, component_types[2].values,
        bottom=component_types[1].get(1, 0), label='Type 2')
plt.bar(component_types[3].index, component_types[3].values, label='Type 3')
plt.bar(component_types[4].index, component_types[4].values,
        bottom=[component_types[3].get(index, 0) for index in component_types[4].index], label='Type 4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.legend()
plt.savefig('out/bar_connectnum-OGnum_type_dist1-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/bar_connectnum-OGnum_type_dist1-2.png')
plt.close()

plt.bar(component_types[3].index, component_types[3].values, label='Type 3', color='C3')
plt.bar(component_types[4].index, component_types[4].values,
        bottom=[component_types[3].get(index, 0) for index in component_types[4].index], label='Type 4', color='C4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.legend()
plt.savefig('out/bar_connectnum-OGnum_type_dist2-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/bar_connectnum-OGnum_type_dist2-2.png')
plt.close()

counts = [component_type.sum() for component_type in component_types]
labels = [f'Type {i}\n({count:,})' for i, count in zip(range(len(component_types)), counts)]
plt.pie(counts, labels=labels, labeldistance=1.2, textprops={'ha': 'center'})
plt.title('Connected components by type')
plt.savefig('out/pie_component_type.png')
plt.close()

"""
DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
./cluster3_graph.py
    ./out/clusters.tsv
"""