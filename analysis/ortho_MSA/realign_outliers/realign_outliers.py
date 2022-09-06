"""Detect alignments with misaligned regions."""

import os
import re
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import skbio
from src.draw import plot_msa_data
from src.utils import read_fasta


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


a = 1E-3  # Coefficient of outlier curve
max_gaps = 1  # Maximum number of gaps in region
min_length = 30  # Minimum length of region
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

records = []
OGids = [path.removesuffix('.afa') for path in os.listdir('../realign_hmmer/out/mafft/') if path.endswith('.afa')]
for OGid in OGids:
    msa = []
    for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa'):
        ppid = re.search(ppid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        msa.append({'ppid': ppid, 'spid': spid, 'seq': seq.upper()})
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])

    ppids = set([record['ppid'] for record in msa])
    spids = set([record['spid'] for record in msa])
    spid2ppids = {spid: [] for spid in spids}
    for record in msa:
        ppid, spid = record['ppid'], record['spid']
        spid2ppids[spid].append(ppid)

    tree = tree_template.shear(spids)
    spid_weights = get_weights(tree)
    ppid_weights = {}
    for spid, ppids in spid2ppids.items():
        weight = spid_weights[spid] / len(ppids)  # Species weight is distributed evenly over all associated genes
        for ppid in ppids:
            ppid_weights[ppid] = weight

    # Count gaps
    gaps = []
    for j in range(len(msa[0]['seq'])):
        gap = sum([1 if msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(msa))])
        gaps.append(gap)

    # Threshold, merge, and size filter to get regions
    binary = ndimage.binary_closing(np.array(gaps) < max_gaps, structure=[1, 1, 1])
    regions = [region for region, in ndimage.find_objects(ndimage.label(binary)[0]) if (region.stop-region.start) >= min_length]

    # Calculate total scores for each sequence over all regions
    scores = {record['ppid']: 0 for record in msa}
    for region in regions:
        for i in range(region.start, region.stop):
            # Build model
            counts = {}
            for record in msa:
                ppid, seq = record['ppid'], record['seq']
                weight = ppid_weights[ppid]
                sym = '-' if seq[i] == '.' else seq[i]  # Convert . to - for counting gaps
                counts[sym] = counts.get(sym, 0) + weight
            total = sum(counts.values())
            model = {sym: np.log(count/total) for sym, count in counts.items()}  # Re-normalize and convert to log space

            # Apply model
            for record in msa:
                ppid, seq = record['ppid'], record['seq']
                sym = '-' if seq[i] == '.' else seq[i]  # Convert . to - for counting gaps
                scores[ppid] += model[sym]

    # Record statistics for each sequence
    values = list(scores.values())
    mean = np.mean(values)
    std = np.std(values)
    iqr = np.quantile(values, 0.75) - np.quantile(values, 0.25)

    for header, score in scores.items():
        x = score - mean
        records.append((x, std, iqr, OGid, msa, regions))

# Plot outputs
if not os.path.exists('out/'):
    os.mkdir('out/')

plt.scatter([record[0] for record in records], [record[1] for record in records], s=10, alpha=0.25, edgecolors='none')
plt.xlabel('(score - alignment mean) of sequence')
plt.ylabel('Standard deviation of scores in alignment')
plt.savefig('out/scatter_std-score.png')
plt.close()

plt.scatter([record[0] for record in records], [record[2] for record in records], s=10, alpha=0.25, edgecolors='none')
plt.xlabel('(score - alignment mean) of sequence')
plt.ylabel('IQR of scores in alignment')
plt.savefig('out/scatter_iqr-score1.png')

ymax = max([record[2] for record in records])
xmin = -(ymax/a)**0.5
xs = np.linspace(xmin, 0, 100)
ys = a*xs**2
plt.plot(xs, ys, color='C1')
plt.savefig('out/scatter_iqr-score2.png')
plt.close()

# Plot unique MSAs with largest deviations
OGids = set()
outliers = sorted([record for record in records if record[0] < -1 and record[2] < a*record[0]**2])  # Use -1 to exclude near-zero floating point rounding errors
for record in outliers:
    # Unpack variables
    OGid, msa, regions = record[3], record[4], record[5]
    if OGid in OGids:
        continue
    OGids.add(OGid)

    # Plot MSA with regions
    line = np.zeros(len(msa[0]['seq']))
    for region in regions:
        line[region] = 1
    plot_msa_data([record['seq'] for record in msa], line, figsize=(15, 6),
                  height_ratio=0.5, hspace=0.2, data_max=1.1, data_min=-0.1, data_labels=['region'],
                  msa_legend=True, legend_kwargs={'bbox_to_anchor': (0.925, 0.5), 'loc': 'center left', 'fontsize': 8, 'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1})
    plt.subplots_adjust(left=0.05, bottom=0.01, right=0.925, top=0.99)
    plt.savefig(f'out/{(len(OGids)-1):03}_{OGid}.png', bbox_inches='tight', dpi=500)
    plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/mafft/*.afa
"""