"""Infer ancestral indel states of IDRs."""

import os
import re
from subprocess import run

import scipy.ndimage as ndimage
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

        # Map deletions to characters
        characters = {}
        for ppid, spid, seq in msa:
            if ppid in ppids:
                binary = [1 if sym in ['-', '.'] else 0 for sym in seq[start:stop]]
                slices = ndimage.find_objects(ndimage.label(binary)[0])
                for s, in slices:
                    try:
                        characters[(s.start, s.stop)].add((ppid, spid))
                    except KeyError:
                        characters[(s.start, s.stop)] = {(ppid, spid)}
        characters = sorted(characters.items(), key=lambda x: (x[0][0], -x[0][1]))  # Fix order of characters

        # Skip region if no indels
        if not characters:
            continue

        # Make character alignment and table
        mca = {(ppid, spid): [] for ppid, spid, _ in msa if ppid in ppids}
        for character, keys1 in characters:
            for key in keys1:
                mca[key].append('1')
            keys0 = set(mca) - keys1
            for key in keys0:
                mca[key].append('0')
        with open(f'out/{prefix}.mfa', 'w') as file:
            for (ppid, spid), charseq in sorted(mca.items(), key=lambda x: x[0][1]):
                seqstring = '\n'.join([''.join(charseq[i:i+80]) for i in range(0, len(charseq), 80)]) + '\n'
                file.write(f'>{spid} {OGid}_{start}-{stop}|{ppid}\n' + seqstring)
        with open(f'out/{prefix}.tsv', 'w') as file:
            file.write('index\tstart\tstop\n')
            for idx, ((s1, s2), _) in enumerate(characters):
                file.write(f'{idx+1}\t{s1}\t{s2}\n')

        # Prune missing species from tree
        spids = set([spid for ppid, spid, _ in msa if ppid in ppids])
        tree = tree_template.shear(spids)
        skbio.io.write(tree, format='newick', into=f'out/{prefix}.nwk')

        run(f'../../../bin/iqtree -s out/{prefix}.mfa -m GTR2+FO+G+ATR -te out/{prefix}.nwk -blfix -pre out/{prefix}', shell=True, check=True)

"""
DEPENDENCIES
../../brownian2/insertion_trim/insertion_trim.py
    ../../brownian2/insertion_trim/out/*.mfa
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
"""