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

# Parse segments
OGid2segments = {}
with open('../../brownian2/aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        if disorder != 'True':
            continue

        record = (int(start), int(stop), set(ppids.split(',')))
        try:
            OGid2segments[OGid].append(record)
        except KeyError:
            OGid2segments[OGid] = [record]

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid, records in OGid2segments.items():
    msa = load_msa(f'../../brownian2/insertion_trim/out/{OGid}.mfa')
    msa = [(re.search(ppid_regex, header).group(1), re.search(spid_regex, header).group(1), seq) for header, seq in msa]

    for start, stop, ppids in records:
        prefix = f'{OGid}_{start}-{stop}'

        # Write region as MSA
        with open(f'out/{prefix}.mfa', 'w') as file:
            for ppid, spid, seq in msa:
                if ppid in ppids:
                    seqstring = '\n'.join([seq[start:stop][i:i+80] for i in range(0, len(seq), 80)]) + '\n'
                    file.write(f'>{spid} {OGid}_{start}-{stop}|{ppid}\n' + seqstring)

        # Prune missing species from tree
        spids = set([spid for ppid, spid, _ in msa if ppid in ppids])
        tree = tree_template.shear(spids)
        skbio.io.write(tree, format='newick', into=f'out/{prefix}.nwk')

        run(f'../../../bin/iqtree -s out/{prefix}.mfa -m 50red_D.txt+I+G -te out/{prefix}.nwk -asr -wslr -pre out/{prefix}', shell=True, check=True)

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
"""