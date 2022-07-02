"""Plot example alignments segmented via posterior decoding."""

import json
import re

import homomorph
import matplotlib.pyplot as plt
import src.draw as draw
from src.utils import read_fasta


class msaBernoulli:
    def __init__(self, ps):
        self.ps = ps

    def pmf(self, x):
        p = self.ps[x[0]]
        if x[1] == 1:
            return p
        else:
            return 1 - p

    def rvs(self, random_state=None):
        # Required for HMM since it has a simulate method
        # Simulations aren't used here, so it's an empty method
        pass


ppid_regex = r'ppid=([A-Za-z0-9_]+)'

# Load model parameters
with open('out/model.json') as file:
    params = json.load(file)


# Load OGids
records = set()
with open('labels.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        records.add((fields['OGid'], fields['ppid']))

# Plot alignments
for OGid, ppid in records:
    # Load MSA
    msa = read_fasta(f'../insertion_trim/out/{OGid}.afa')
    seq = [seq for header, seq in msa if re.search(ppid_regex, header).group(1) == ppid][0]

    # Create Bernoulli sequence
    ps = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        p = sum(col) / len(col)
        ps.append(p)

    # Create emission sequence
    emits = []
    for j, sym in enumerate(seq):
        if sym in ['-', '.']:
            emits.append((j, 1))
        else:
            emits.append((j, 0))

    # Instantiate model
    e_dists_rv = {'0': msaBernoulli(ps), '1': msaBernoulli([params['e_param'] for _ in range(len(ps))])}
    model = homomorph.HMM(params['t_dists'], e_dists_rv, params['start_dist'])

    # Decode states and plot
    fbs = model.forward_backward(emits)
    draw.plot_msa_data([seq for _, seq in msa], fbs['1'],
                       msa_labels=[ppid if ppid in header else '' for header, _ in msa], msa_labelsize=4,
                       figsize=(15, 6))
    plt.savefig(f'out/{OGid}_{ppid}_wide.png', bbox_inches='tight')
    plt.close()

    draw.plot_msa_data([seq for _, seq in msa], fbs['1'],
                       msa_labels=[ppid if ppid in header else '' for header, _ in msa], msa_labelsize=4,
                       figsize=(8, 8))
    plt.savefig(f'out/{OGid}_{ppid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../insertion_trim/extract.py
    ../trim_extract/out/*.afa
./fit.py
    ./out/model.json
./labels.tsv
"""