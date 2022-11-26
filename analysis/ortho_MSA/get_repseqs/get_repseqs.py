"""Select representative sequences from alignments."""

import os
import re

import numpy as np
import skbio
from src.utils import get_brownian_weights, read_fasta


def get_score(msa_model, seq):
    score = 0
    for column_model, sym in zip(msa_model, seq):
        sym = 1 if sym != '-' else 0
        score += column_model[sym]
    return score


gnid_regex = r'gnid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'
pseudocount = 0.01

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path.removesuffix('.afa') for path in os.listdir('../align_fastas/out/') if path.endswith('.afa')]
for OGid in OGids:
    input_msa = []
    for header, seq in read_fasta(f'../align_fastas/out/{OGid}.afa'):
        gnid = re.search(gnid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        input_msa.append({'header': header, 'gnid': gnid, 'spid': spid, 'seq': seq})

    gnids = set([record['gnid'] for record in input_msa])
    spids = set([record['spid'] for record in input_msa])
    gnid2records = {gnid: [] for gnid in gnids}
    spid2gnids = {spid: set() for spid in spids}
    for record in input_msa:
        gnid, spid = record['gnid'], record['spid']
        gnid2records[gnid].append(record)
        spid2gnids[spid].add(gnid)

    tree = tree_template.shear(spids)
    spid_weights = {tip.name: weight for tip, weight in get_brownian_weights(tree)}
    gnid_weights = {}
    for spid, gnids in spid2gnids.items():
        weight = spid_weights[spid] / len(gnids)  # Species weight is distributed evenly over all associated genes
        for gnid in gnids:
            gnid_weights[gnid] = weight

    msa_model = []
    for j in range(len(input_msa[0]['seq'])):
        counts = {0: pseudocount/2, 1: pseudocount/2}
        for gnid, records in gnid2records.items():
            weight = gnid_weights[gnid]
            any_aligned = any(record['seq'][j] != '-' for record in records)
            if any_aligned:
                counts[1] += weight
            else:
                counts[0] += weight
        total = sum(counts.values())
        column_model = {sym: np.log(count / total) for sym, count in counts.items()}  # Re-normalize and convert to log space
        msa_model.append(column_model)

    output_msa = []
    for records in gnid2records.values():
        record = max(records, key=lambda x: get_score(msa_model, x['seq']))
        output_msa.append(record)
    output_msa = sorted(output_msa, key=lambda x: x['spid'])

    # Remove excess gaps
    slices, idx = [], None
    for j in range(len(output_msa[0]['seq'])):
        for i in range(len(output_msa)):
            sym = output_msa[i]['seq'][j]
            if sym not in ['-', '.']:
                if idx is None:  # Store position only if new slice is not started
                    idx = j
                break
        else:
            if idx is not None:
                slices.append(slice(idx, j))
                idx = None
    if idx is not None:  # Add final slice to end
        slices.append(slice(idx, len(output_msa[0]['seq'])))

    with open(f'out/{OGid}.afa', 'w') as file:
        for record in output_msa:
            header, seq1 = record['header'], record['seq']
            seq2 = ''.join([seq1[s] for s in slices])  # Remove gap only columns
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')

"""
NOTES
Since the model for each gene is constructed from a kind of union of its sequences, it's possible for some positions to
have a 100% fraction aligned even if some sequences are unaligned. Specifically, this will occur if every gene has at
least one sequence aligned at that positions and at least one gene has a sequence that is not aligned at that position.
To prevent these sequencing from scoring -infinity, the counts are initialized with pseudocounts. The counts at each
position are already weighted to sum to 1, so the pseudocounts are read as a ratio of the weight of the prior to the
weight of the data. (The pseudocount is split evenly across gap and non-gap outcomes to not double the weight of the
prior.)

DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../align_fastas/align_fastas.py
    ../align_fastas/out/*.afa
"""