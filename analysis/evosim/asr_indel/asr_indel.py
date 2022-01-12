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

        # Make list of gaps for each sequence
        ids2characters = {}
        for ppid, spid, seq in msa:
            if ppid in ppids:
                binary = [1 if sym in ['-', '.'] else 0 for sym in seq[start:stop]]
                slices = ndimage.find_objects(ndimage.label(binary)[0])
                ids2characters[(ppid, spid)] = [(s.start, s.stop) for s, in slices]
        character_set = set().union(*ids2characters.values())
        character_set = sorted(character_set, key=lambda x: (x[0], -x[1]))  # Fix order of characters

        # Skip region if no indels
        if not character_set:
            continue

        # Make character alignment and table
        mca = {}
        for ids, characters in ids2characters.items():
            charseq = []
            for (start1, stop1) in character_set:
                for (start2, stop2) in characters:
                    if (start2 <= start1) and (stop2 >= stop1):
                        charseq.append('1')
                        break
                else:
                    charseq.append('0')
            mca[ids] = charseq
        with open(f'out/{prefix}.mfa', 'w') as file:
            for (ppid, spid), charseq in sorted(mca.items(), key=lambda x: x[0][1]):
                seqstring = '\n'.join([''.join(charseq[i:i+80]) for i in range(0, len(charseq), 80)]) + '\n'
                file.write(f'>{spid} {OGid}_{start}-{stop}|{ppid}\n' + seqstring)
        with open(f'out/{prefix}.tsv', 'w') as file:
            file.write('index\tstart\tstop\n')
            for idx, (s1, s2) in enumerate(character_set):
                file.write(f'{idx}\t{s1}\t{s2}\n')

        # Prune missing species from tree
        spids = set([spid for ppid, spid, _ in msa if ppid in ppids])
        tree = tree_template.shear(spids)
        skbio.io.write(tree, format='newick', into=f'out/{prefix}.nwk')

        run(f'../../../bin/iqtree -s out/{prefix}.mfa -m GTR2+FO+G+ASC -te out/{prefix}.nwk -blfix -keep-ident -pre out/{prefix}', shell=True, check=True)

"""
NOTES
Test reconstructions showed that treating every character as independent underestimates the probability of gaps. I
believe the issue is if a region with consistent gaps has "ragged" ends, each gap with a unique start or stop position
is coded as a separate character. In the most extreme case, a gap common to every sequence but one may differ by at
least one at each stop position, like a staircase. Thus, the nested structure of the gaps is not reflected in the
character codings.

DEPENDENCIES
../../brownian2/insertion_trim/insertion_trim.py
    ../../brownian2/insertion_trim/out/*.mfa
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
"""