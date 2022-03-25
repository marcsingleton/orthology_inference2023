"""Fit rates to GO term logistic regression models."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

pdidx = pd.IndexSlice

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load regions
rows = []
with open('../../brownian2/aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        ppid = re.search(r'(FBpp[0-9]+)', ppids).group(1)
        rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'disorder': disorder == 'True', 'gnid': ppid2gnid[ppid]})
regions = pd.DataFrame(rows)

# Load GOids
GOids, gnid2GOids = set(), {}
with open('../filter_GAF/out/GAF_drop.tsv') as file:
    header = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(header, line.rstrip('\n').split('\t'))}
        gnid, GOid = fields['gnid'], fields['GOid']
        GOids.add(GOid)
        try:
            gnid2GOids[gnid].add(GOid)
        except KeyError:
            gnid2GOids[gnid] = {GOid}

contrasts = pd.read_table('../../brownian2/get_contrasts/out/contrasts.tsv')
df1 = regions.merge(contrasts, how='right', on=['OGid', 'start', 'stop']).set_index(['OGid', 'start', 'stop', 'disorder', 'gnid', 'contrast_id'])

rates = ((df1**2).groupby(['OGid', 'start', 'stop', 'disorder', 'gnid']).mean())
disorder = rates.loc[pdidx[:, :, :, True], :]
order = rates.loc[pdidx[:, :, :, False], :]

if not os.path.exists('out/'):
    os.mkdir('out/')

rows = []
for data, label, color in [(disorder, 'disorder', 'C0'), (order, 'order', 'C1'), (rates, 'combined', 'C2')]:
    corrs = data.corr()
    plt.figure(figsize=(8, 6))
    plt.imshow(corrs, vmin=0, vmax=1)
    plt.xticks(range(len(corrs.index)), corrs.index, fontsize=6, rotation=60, rotation_mode='anchor',
               horizontalalignment='right', verticalalignment='center')
    plt.yticks(range(len(corrs.index)), corrs.index, fontsize=6)
    plt.colorbar(shrink=0.75)
    plt.savefig(f'out/heatmap_corr_{label}.png', bbox_inches='tight')
    plt.close()

    pca = PCA(n_components=10)
    transform = pca.fit_transform(data.to_numpy())[:, :5]

    plt.bar(range(1, len(pca.explained_variance_ratio_) + 1), pca.explained_variance_ratio_, label=label, color=color)
    plt.xlabel('Principal component')
    plt.ylabel('Explained variance ratio')
    plt.legend()
    plt.savefig(f'out/bar_scree_{label}.png')
    plt.close()

    for i, GOid in enumerate(GOids):
        y_true = [GOid in gnid2GOids.get(gnid, set()) for gnid in data.index.get_level_values('gnid')]
        w = (len(y_true) - sum(y_true)) / sum(y_true)
        weights = [w if y else 1 for y in y_true]
        logit = LogisticRegression(max_iter=500, penalty='none')
        logit.fit(transform, y_true, sample_weight=weights)

        y_pred = logit.predict(transform)
        acc = sum([y1 == y2 for y1, y2 in zip(y_true, y_pred)]) / len(y_true)
        sens = sum([(y1 == y2) and y1 for y1, y2 in zip(y_true, y_pred)]) / sum(y_true)
        spec = sum([(y1 == y2) and not y1 for y1, y2 in zip(y_true, y_pred)]) / (len(y_true) - sum(y_true))
        prec = sum([(y1 == y2) and y1 for y1, y2 in zip(y_true, y_pred)]) / sum(y_pred)

        rows.append({'GOid': GOid, 'label': label,
                     'accuracy': acc, 'sensitivty': sens, 'specificity': spec, 'precision': prec,
                     **{f'beta{i}': beta for i, beta in enumerate(logit.coef_[0])}})
results = pd.DataFrame(rows)
results.to_csv('out/models.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../brownian2/aucpred_filter/aucpred_filter.py
    ../../brownian2/aucpred_filter/out/regions_30.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../filter_GAF/filter_GAF.py
    ../filter_GAF/out/GAF_drop.tsv
"""