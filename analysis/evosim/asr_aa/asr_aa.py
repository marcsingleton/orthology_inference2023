"""Infer ancestral amino acid distributions of IDRs."""

import os
import re
from subprocess import run

import skbio
from src.utils import read_fasta

ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)

OGids = set()
with open('../../brownian2/aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.rstrip('\n').split('\t')
        OGids.add(OGid)

OGid2regions = {}
with open('../../brownian2/aucpred_regions/out/regions.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.rstrip('\n').split('\t')
        try:
            OGid2regions[OGid].append((int(start), int(stop), True if disorder == 'True' else False))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), True if disorder == 'True' else False)]

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid in OGids:
    msa = read_fasta(f'../../brownian2/insertion_trim/out/{OGid}.afa')
    msa = [(re.search(ppid_regex, header).group(1), re.search(spid_regex, header).group(1), seq) for header, seq in msa]

    # Check regions and merge if necessary
    regions = OGid2regions[OGid]
    disorder_length = sum([stop-start for start, stop, disorder in regions if disorder])
    order_length = sum([stop-start for start, stop, disorder in regions if not disorder])
    if disorder_length >= 30 and order_length >= 30:
        disorder_regions = [f'{start+1}-{stop}' for start, stop, disorder in regions if disorder]
        order_regions = [f'{start+1}-{stop}' for start, stop, disorder in regions if not disorder]
    elif disorder_length >= 30:
        disorder_regions = [f'1-{len(msa[0][2])}']
        order_regions = []
    elif order_length >= 30:
        disorder_regions = []
        order_regions = [f'1-{len(msa[0][2])}']
    else:
        continue

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
    with open(f'out/{OGid}.afa', 'w') as file:
        for ppid, spid, seq in msa:
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            file.write(f'>{spid} {ppid}\n{seqstring}\n')

    # Make NEXUS partition file
    with open(f'out/{OGid}.nex', 'w') as file:
        partitions = []
        file.write('#nexus\nbegin sets;\n')
        if disorder_regions:
            disorder_string = ' '.join(disorder_regions)
            partitions.append('../config/50R_disorder.paml+I+G:disorder')
            file.write(f'    charset disorder = {disorder_string};\n')
        if order_regions:
            order_string = ' '.join(order_regions)
            partitions.append('LG+I+G:order')
            file.write(f'    charset order = {order_string};\n')
        partition_string = ', '.join(partitions)
        file.write(f'    charpartition mine = {partition_string};\nend;\n')

    # Prune missing species from tree
    spids = {spid for _, spid, _ in msa}
    tree = tree_template.shear(spids)
    skbio.io.write(tree, format='newick', into=f'out/{OGid}.nwk')

    run(f'../../../bin/iqtree -s out/{OGid}.afa -spp out/{OGid}.nex -te out/{OGid}.nwk -keep-ident -pre out/{OGid}', shell=True, check=True)

"""
NOTES
The documentation for tree and branch arguments to IQ-TREE does not entirely describe what is optimized and what is
fixed. For example, using -te fixes the tree topology, but it scales the branch lengths individually. Using -t in
combination with -blscale will also fix the tree topology but it will scale all branch lengths by a constant value.
Finally, using -t alone will perform a full tree search and branch length optimization starting at the given tree (as
described in the documentation.)

When the number of parameters exceeds the number of columns, IQ-TREE cautions that the parameter estimates are
unreliable. Since a majority of the parameters are branch lengths, I attempted to reduce the model complexity by using
the -blscale option. Unfortunately, the -blscale option does not play nicely with partitioned models, and IQ-TREE
crashes when both are used together. While this is not ideal, it is likely of minor importance. Since the overall goal
is to generate a sample of plausible ancestral sequences, the exact parameter estimates are unimportant. Furthermore,
the model itself (with its heuristic treatment of indels) is likely a poor approximation for the process that generated
the sequences. I instead consider this model as a "better" consensus that incorporates evolutionary relationships and
prior expectations about the frequencies and exchangeabilities of amino acids.

To prevent poor fits from a lack of data, a partition is only created if there are a minimum of 30 columns in that 
class. If one of the classes has 30 columns but the other does not, the small class is consolidated into the large one.
If neither class has 30 columns, the alignment is skipped. These rules ensure that any alignment with region in the 
final set of regions is reconstructed. They also ensure any regions in the final set are reconstructed with a model that
is (largely) fit to its class.

DEPENDENCIES
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
../../brownian2/aucpred_regions/get_regions.py
    ../../brownian2/aucpred_regions/out/regions.tsv
../../brownian2/insertion_trim/extract.py
    ../../brownian2/insertion_trim/out/*.afa
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../config/50R_disorder.paml
"""