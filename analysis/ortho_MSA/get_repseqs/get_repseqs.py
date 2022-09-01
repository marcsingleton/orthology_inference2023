"""Select representative sequences from alignments."""

import os
import re
from itertools import product

import numpy as np
import skbio
from src.utils import read_fasta


def get_score(msa_model, seq):
    score = 0
    for column_model, sym in zip(msa_model, seq):
        sym = 1 if sym != '-' else 0
        score += column_model[sym]
    return score


def get_weights(tree):
    spids = [tip.name for tip in tree.tips()]
    spid2idx = {spid: i for i, spid in enumerate(spids)}
    idx2spid = {i: spid for i, spid in enumerate(spids)}

    # Accumulate tip names up to root
    for node in tree.postorder():
        if node.is_tip():
            node.spids = {node.name}
        else:
            node.spids = set.union(*[child.spids for child in node.children])

    # Fill in covariance matrix from root to tips
    tree.root_length = 0
    cov = np.zeros((len(spids), len(spids)))
    for node in tree.traverse(include_self=True):
        for child in node.children:
            child.root_length = node.root_length + child.length
        if not node.is_tip():
            child1, child2 = node.children
            idxs1, idxs2 = [spid2idx[spid] for spid in child1.spids], [spid2idx[spid] for spid in child2.spids]
            for idx1, idx2 in product(idxs1, idxs2):
                cov[idx1, idx2] = node.root_length
                cov[idx2, idx1] = node.root_length
        else:
            idx = spid2idx[node.name]
            cov[idx, idx] = node.root_length

    # Compute weights
    # The formula below is from the appendix of J. Mol. Biol. (1989) 207. 647-653.
    # It assumes continuous traits evolve by a Brownian motion model
    # The MLE for the root value is then a weighted average of the observed tip values with the weights given below
    inv = np.linalg.inv(cov)
    row_sum = inv.sum(axis=1)
    total_sum = inv.sum()
    ws = row_sum / total_sum
    return {idx2spid[idx]: w for idx, w in enumerate(ws)}


gnid_regex = r'gnid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'
pseudocount = 0.01

tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path.removesuffix('.afa') for path in os.listdir('../align_fastas/out/') if path.endswith('.afa')]
for OGid in OGids:
    input_msa = []
    for header, seq in read_fasta(f'../align_fastas/out/{OGid}.afa'):
        gnid = re.search(gnid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        input_msa.append({'header': header, 'gnid': gnid, 'spid': spid, 'seq': seq})
    input_msa = sorted(input_msa, key=lambda x: tip_order[x['spid']])

    gnids = set([record['gnid'] for record in input_msa])
    spids = set([record['spid'] for record in input_msa])
    gnid2records = {gnid: [] for gnid in gnids}
    spid2gnids = {spid: [] for spid in spids}
    for record in input_msa:
        gnid, spid = record['gnid'], record['spid']
        gnid2records[gnid].append(record)
        spid2gnids[spid].append(gnid)

    tree = tree_template.shear(spids)
    spid_weights = get_weights(tree)
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
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../align_fastas/align_fastas.py
    ../align_fastas/out/*.afa
"""