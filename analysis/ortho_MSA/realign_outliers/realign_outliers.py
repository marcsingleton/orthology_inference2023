"""Detect alignments with misaligned regions."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import skbio
from src.draw import plot_msa_data
from src.utils import read_fasta, get_brownian_weights

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

a = 2E-3  # Coefficient of outlier curve
b = -1  # Intercept of outlier curve
max_gaps = 1  # Maximum number of gaps in region
min_length = 30  # Minimum length of region

tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

records = []
OGids = [path.removesuffix('.afa') for path in os.listdir('../realign_fastas/out/') if path.endswith('.afa')]
for OGid in OGids:
    # Load MSA
    msa = []
    for header, seq in read_fasta(f'../realign_fastas/out/{OGid}.afa'):
        ppid = re.search(ppid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        msa.append({'ppid': ppid, 'spid': spid, 'seq': seq.upper()})
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])

    # Calculate weights
    spid_set = set([record['spid'] for record in msa])
    spid2ppids = {spid: [] for spid in spid_set}
    for record in msa:
        ppid, spid = record['ppid'], record['spid']
        spid2ppids[spid].append(ppid)

    tree = tree_template.shear(spid_set)
    spid_weights = {tip.name: weight for tip, weight in get_brownian_weights(tree)}
    ppid_weights = {}
    for spid, ppids in spid2ppids.items():
        weight = spid_weights[spid] / len(ppids)  # Species weight is distributed evenly over all associated proteins
        for ppid in ppids:
            ppid_weights[ppid] = weight

    # Count gaps
    gaps = []
    for j in range(len(msa[0]['seq'])):
        gap = sum([1 if msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(msa))])
        gaps.append(gap)

    # Threshold, merge, and size filter to get regions
    binary = ndimage.binary_closing(np.array(gaps) < max_gaps, structure=[1, 1, 1])
    slices = [s for s, in ndimage.find_objects(ndimage.label(binary)[0]) if (s.stop-s.start) >= min_length]

    # Calculate total scores for each sequence over all regions
    scores = {record['ppid']: 0 for record in msa}
    for s in slices:
        for i in range(s.start, s.stop):
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
        records.append((x, std, iqr, OGid, msa, slices))

# Plot outputs
if not os.path.exists('out/'):
    os.mkdir('out/')

plt.hexbin([record[0] for record in records], [record[1] for record in records], bins='log', gridsize=75, mincnt=1, linewidth=0)
plt.xlabel('(score - alignment mean) of sequence')
plt.ylabel('Standard deviation of scores in alignment')
plt.colorbar()
plt.savefig('out/hexbin_std-score.png')
plt.close()

plt.hexbin([record[0] for record in records], [record[2] for record in records], bins='log', gridsize=75, mincnt=1, linewidth=0)
plt.xlabel('(score - alignment mean) of sequence')
plt.ylabel('IQR of scores in alignment')
plt.colorbar()
plt.savefig('out/hexbin_iqr-score1.png')

ymax = max([record[2] for record in records])
xmin = -(ymax/a)**0.5
xs = np.linspace(xmin, 0, 100)
ys = a*xs**2 + b
plt.plot(xs, ys, color='C1')
plt.savefig('out/hexbin_iqr-score2.png')
plt.close()

# Plot unique MSAs with largest deviations
OGids = set()
outliers = sorted([record for record in records if record[0] < 0 and record[2] < a*record[0]**2 + b])
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
    plot_msa_data([record['seq'] for record in msa], line,
                  figsize=(15, 6),
                  height_ratio=0.5, hspace=0.45, left=0.015, right=0.94, top=0.99, bottom=0.03,
                  data_max=1.05, data_min=-0.05, data_labels=['region'],
                  msa_legend=True, legend_kwargs={'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                  'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1})
    plt.savefig(f'out/{(len(OGids)-1):03}_{OGid}.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../realign_fastas/realign_fastas.py
    ../realign_fastas/out/*.afa
"""