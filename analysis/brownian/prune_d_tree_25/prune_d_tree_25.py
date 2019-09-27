""""Prune unneeded lineages from Drosophila tree and save image of pruned tree."""

import Bio.Phylo as Phylo
import matplotlib.pyplot as plt

remain = set(['ananassae', 'erecta', 'grimshawi', 'melanogaster', 'mojavensis', 'persimilis', 'pseudoobscura', 'virilis', 'yakuba', 'willistoni'])

tree = Phylo.read('../../../data/d_tree_25/drosophila-25spec-tree.tre', 'newick')
leaf_names = set([terminal.name for terminal in tree.get_terminals()])

remove = leaf_names - remain
for leaf in remove:
    tree.prune(leaf)

Phylo.write(tree, 'drosophila-10spec-tree.nwk', 'newick')

fig, ax = plt.subplots(figsize=(7, 4))
Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
plt.savefig('drosophila-10spec-tree.png')
