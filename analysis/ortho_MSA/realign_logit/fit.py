"""Fit logistic regression using labelled segments."""

import json
import multiprocessing as mp
import os
from itertools import product

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import skbio
from src.ortho_MSA.trim import trim_insertions
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, log_loss
from src.utils import read_fasta


def fit_model(OGid2msa,
              con_frac, con_window, con_minlen,
              gap_num, gap_rate, gap_minsig,
              nongap_frac, nongap_minlen,
              gp_sigma, gd_window, indel1_rate, indel2_rate,
              weights, threshold,
              matrix):
    # Calculate regressors for each OG and combine into df
    rows = []
    for OGid, (msa, scores, gaps_array) in OGid2msa.items():
        rows.extend(get_regressors(OGid, msa, scores, gaps_array,
                                   con_frac, con_window, con_minlen,
                                   gap_num, gap_rate, gap_minsig,
                                   nongap_frac, nongap_minlen,
                                   gp_sigma, gd_window, indel1_rate, indel2_rate,
                                   weights, threshold,
                                   matrix))
    df = labels.merge(pd.DataFrame(rows), on=['OGid', 'start', 'stop', 'index'], how='inner')
    df['trimmed3'] = df['trimmed1'] | df['trimmed2']

    # Fit model
    regressors = ['length', 'support', 'gap_propensity', 'gap_diversity', 'indel_bias1', 'indel_bias2', 'terminal']
    X = df[regressors]
    y_true = df['trimmed3'].astype(bool)
    logit = LogisticRegression(max_iter=500, penalty='none')
    logit.fit(X, y_true)

    # Get model parameters
    d = {'gp_sigma': gp_sigma, 'gd_window': gd_window, 'indel1_rate': indel1_rate, 'indel2_rate': indel2_rate}

    d['bias'] = logit.intercept_[0]
    d.update({key: value for key, value in zip(regressors, logit.coef_[0])})

    y_pred = logit.predict(X)
    d['log_loss'] = log_loss(y_true, y_pred)
    d.update({key: value for key, value in zip(['tn', 'fp', 'fn', 'tp'], confusion_matrix(y_true, y_pred).ravel())})

    return d


def get_regressors(OGid, msa, scores, gaps_array,
                   con_frac, con_window, con_minlen,
                   gap_num, gap_rate, gap_minsig,
                   nongap_frac, nongap_minlen,
                   gp_sigma, gd_window, indel1_rate, indel2_rate,
                   weights, threshold,
                   matrix):
    binary = ndimage.binary_closing(scores <= len(msa) * con_frac)
    mask = ndimage.label(binary, structure=con_window * [1])[0]
    regions = [region for region, in ndimage.find_objects(mask) if region.stop - region.start >= con_minlen]
    nterm = min([region.start for region in regions])
    cterm = max([region.stop for region in regions])

    ds = []
    _, trims = trim_insertions(msa, scores, gaps_array,
                               gap_num, gap_rate, gap_minsig,
                               nongap_frac, nongap_minlen,
                               gp_sigma, gd_window, indel1_rate, indel2_rate,
                               weights, threshold,
                               matrix)
    for trim in trims:
        if trim['trimmed'] is not None:
            region = trim['region']
            if region.stop <= nterm or region.start >= cterm:
                terminal = 1
            else:
                terminal = 0

            keys = ['index', 'length', 'support', 'gap_propensity', 'gap_diversity', 'indel_bias1', 'indel_bias2']
            d = {'OGid': OGid, 'start': region.start, 'stop': region.stop, 'terminal': terminal,
                 **{key: trim[key] for key in keys}}
            ds.append(d)
    return ds


# Load parameters
num_processes = int(os.environ['SLURM_NTASKS'])

with open('../config/trim_params.json') as file:
    tp = json.load(file)
matrix = {}
with open('../config/BLOSUM62.txt') as file:
    file.readline()  # Skip header
    syms = file.readline().split()
    for i, line in enumerate(file):
        for j, value in enumerate(line.split()[1:]):
            matrix[(syms[i], syms[j])] = int(value)
weights = {key: 1 for key in ['bias', 'length', 'support', 'gap_propensity', 'gap_diversity', 'indel_bias1', 'indel_bias2']}

labels = pd.read_table('out/segments_label.tsv').dropna().drop('length', axis=1)

gp_range = np.linspace(1, 4, 7)
gd_range = np.linspace(5, 25, 5)
ir1_range = np.linspace(0.01, 0.09, 9)
ir2_range = np.linspace(0.01, 0.09, 9)
ranges = product(gp_range, gd_range, ir1_range, ir2_range)

if __name__ == '__main__':
    # Load MSAs
    OGid2msa = {}
    for OGid in labels['OGid'].drop_duplicates():
        try:
            msa = read_fasta(f'../align_fastas1/out/{OGid}.mfa')
        except FileNotFoundError:
            msa = read_fasta(f'../align_fastas2-2/out/{OGid}.mfa')

        gaps_array = np.full((len(msa), len(msa[0][1])), False)
        for i, (_, seq) in enumerate(msa):
            for j, sym in enumerate(seq):
                if sym == '-':
                    gaps_array[i, j] = True
        scores = gaps_array.sum(axis=0)
        msa = skbio.TabularMSA([skbio.Protein(seq, metadata={'description': header}) for header, seq in msa])
        OGid2msa[OGid] = (msa, scores, gaps_array)

    # Construct and apply arguments
    args_list = []
    for gp_sigma, gd_window, indel1_rate, indel2_rate in ranges:
        args = [OGid2msa,
                tp['con_frac'], tp['con_window'], tp['con_minlen'],
                tp['gap_num'], tp['gap_rate'], tp['gap_minsig'],
                tp['nongap_frac'], tp['nongap_minlen'],
                gp_sigma, gd_window, indel1_rate, indel2_rate,
                weights, tp['threshold'],
                matrix]
        args_list.append(args)
    with mp.Pool(processes=num_processes) as pool:
        rows = pool.starmap(fit_model, args_list)

    # Save results
    if not os.path.exists('out/'):
        os.mkdir('out/')

    df = pd.DataFrame(rows)
    df.to_csv('out/models.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../config/BLOSUM62.txt
../config/trim_params.json
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
./extract.py
"""