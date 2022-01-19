"""Infer ancestral amino acid distributions of IDRs."""

import os
import re
from subprocess import run

import skbio


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'
tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)

OGid2regions = {}
with open('../../brownian2/aucpred_regions/out/regions.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.split()
        try:
            OGid2regions[OGid].append((int(start), int(stop), True if disorder == 'True' else False))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), True if disorder == 'True' else False)]

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path[:-4] for path in os.listdir('../../brownian2/insertion_trim/out/') if path.endswith('.mfa')]
for OGid in OGids:
    msa = load_msa(f'../../brownian2/insertion_trim/out/{OGid}.mfa')
    msa = [(re.search(ppid_regex, header).group(1), re.search(spid_regex, header).group(1), seq) for header, seq in msa]

    # Skip MSA if is invariant
    is_invariant = True
    for j in range(len(msa[0][2])):
        sym0 = msa[0][2][j]
        if any([msa[i][2][j] != sym0 for i in range(1, len(msa))]):
            is_invariant = False
            break
    if is_invariant:
        continue

    # Write region as MSA
    with open(f'out/{OGid}.mfa', 'w') as file:
        for ppid, spid, seq in msa:
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
            file.write(f'>{spid} {ppid}\n' + seqstring)

    # Make NEXUS partition file
    regions = OGid2regions[OGid]
    disorder_regions = [f'{start+1}-{stop}' for start, stop, disorder in regions if disorder]
    order_regions = [f'{start+1}-{stop}' for start, stop, disorder in regions if not disorder]
    with open(f'out/{OGid}.nex', 'w') as file:
        partitions = []
        file.write('#nexus\nbegin sets;\n')
        if disorder_regions:
            disorder_string = ' '.join(disorder_regions)
            partitions.append('50red_D.txt+I+G:disorder')
            file.write(f'    charset disorder = {disorder_string};\n')
        if order_regions:
            order_string = ' '.join(order_regions)
            partitions.append('WAG+I+G:order')
            file.write(f'    charset order = {order_string};\n')
        partition_string = ', '.join(partitions)
        file.write(f'    charpartition mine = {partition_string};\nend;\n')

    # Prune missing species from tree
    spids = set([spid for _, spid, _ in msa])
    tree = tree_template.shear(spids)
    skbio.io.write(tree, format='newick', into=f'out/{OGid}.nwk')

    run(f'../../../bin/iqtree -s out/{OGid}.mfa -spp out/{OGid}.nex -te out/{OGid}.nwk -keep-ident -pre out/{OGid}', shell=True, check=True)

"""
NOTES
The documentation for tree and branch arguments to IQ-TREE does not entirely describe what is optimized and what is
fixed. For example, using -te fixes the tree topology, but it scales the branch lengths individually. Using -t in
combination with -blscale will also fix the tree topology but it will scale all branch lengths by a constant value.
Finally, using -t alone will perform a full tree search and branch length optimization starting at the given tree (as
described in the documentation.)

DEPENDENCIES
../../brownian2/insertion_trim/insertion_trim.py
    ../../brownian2/insertion_trim/out/*.mfa
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
./50red_D.txt
"""