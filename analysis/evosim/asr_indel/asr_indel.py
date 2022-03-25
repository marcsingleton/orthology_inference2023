"""Infer ancestral indel states of IDRs."""

import os
import re
from subprocess import run

import scipy.ndimage as ndimage
import skbio
from src.utils import read_fasta


def is_nested(character, characters):
    """Return if character is nested in one of intervals in characters."""
    start1, stop1 = character
    for start2, stop2 in characters:
        if (start2 <= start1) and (stop2 >= stop1):
            return True
    return False


ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)

OGids = set()
with open('../../brownian2/aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        OGids.add(OGid)

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

for OGid in OGids:
    msa = read_fasta(f'../../brownian2/insertion_trim/out/{OGid}.mfa')
    msa = [(re.search(ppid_regex, header).group(1), re.search(spid_regex, header).group(1), seq) for header, seq in msa]

    # Check regions (continuing only if alignment is fit by asr_aa.py)
    regions = OGid2regions[OGid]
    disorder_length = sum([stop-start for start, stop, disorder in regions if disorder])
    order_length = sum([stop-start for start, stop, disorder in regions if not disorder])
    if disorder_length <= 30 and order_length <= 30:
        continue

    # Make list of gaps for each sequence
    ids2characters = {}
    for ppid, spid, seq in msa:
        binary = [1 if sym in ['-', '.'] else 0 for sym in seq]
        slices = ndimage.find_objects(ndimage.label(binary)[0])
        ids2characters[(ppid, spid)] = [(s.start, s.stop) for s, in slices]
    character_set = set().union(*ids2characters.values())
    character_set = sorted(character_set, key=lambda x: (x[0], -x[1]))  # Fix order of characters

    # Skip region if no indels
    if not character_set:
        continue

    # Make character alignment
    mca = []
    for (ppid, spid), characters in ids2characters.items():
        charseq = []
        for character in character_set:
            if is_nested(character, characters):
                charseq.append('1')
            else:
                charseq.append('0')
        mca.append((ppid, spid, charseq))
    mca = sorted(mca, key=lambda x: x[1])

    # Identify invariant characters
    is_invariants = []
    for j in range(len(mca[0][2])):
        is_invariant = True
        for i in range(len(mca)):
            if mca[i][2][j] == '0':
                is_invariant = False
                break
        is_invariants.append(is_invariant)

    # Write character table to file
    idx = 0
    with open(f'out/{OGid}.tsv', 'w') as file:
        file.write('index\tstart\tstop\n')
        for is_invariant, (start, stop) in zip(is_invariants, character_set):
            if is_invariant:
                file.write(f'-1\t{start}\t{stop}\n')
            else:
                file.write(f'{idx}\t{start}\t{stop}\n')
                idx += 1

    # Skip model fit if all characters are invariant
    if all(is_invariants):
        continue

    # Write alignment to file
    with open(f'out/{OGid}.mfa', 'w') as file:
        for ppid, spid, charseq in mca:
            charseq = [sym for is_invariant, sym in zip(is_invariants, charseq) if not is_invariant]  # Filter invariant characters
            seqstring = '\n'.join([''.join(charseq[i:i+80]) for i in range(0, len(charseq), 80)]) + '\n'
            file.write(f'>{spid} {OGid}_{start}-{stop}|{ppid}\n' + seqstring)

    # Prune missing species from tree
    spids = {spid for _, spid, _ in msa}
    tree = tree_template.shear(spids)
    skbio.io.write(tree, format='newick', into=f'out/{OGid}.nwk')

    run(f'../../../bin/iqtree -s out/{OGid}.mfa -m GTR2+FO+G+ASC -te out/{OGid}.nwk -keep-ident -pre out/{OGid}', shell=True, check=True)

"""
NOTES
Test reconstructions showed that treating every character as independent underestimates the probability of gaps. I
believe the issue is if a region with consistent gaps has "ragged" ends, each gap with a unique start or stop position
is coded as a separate character. In the most extreme case, a gap common to every sequence but one may differ by at
least one at each stop position, like a staircase. Thus, the nested structure of the gaps is not reflected in the
character codings.

DEPENDENCIES
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
../../brownian2/aucpred_regions/get_regions.py
    ../../brownian2/aucpred_regions/out/regions.tsv
../../brownian2/insertion_trim/extract.py
    ../../brownian2/insertion_trim/out/*.mfa
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
"""