"""Segment alignments by IUPRED2a and combine by union."""

import gzip
import os
import numpy as np
import pandas as pd
import subprocess
import tempfile
from Bio import AlignIO
from scipy import ndimage

# Input variables
path = '../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv'
dir = '../../EggNOGv5_validation/filter_unknown_realign/out/align/'
thresh = 0.5

block_data = []  # Block data with raw sequences
block_num = 0  # Counter for numbering blocks
seg_num = 0  # Counter for numbering rows

with open(path) as file:
    for line in file:
        fields = line.rstrip('\n').split('\t')
        ali_id = fields[1]
        with gzip.open(dir + ali_id + '.raw_alg.faa.gz', 'rt') as file_MSA:  # read, text mode
            MSA = AlignIO.read(file_MSA, 'fasta')
            MSA_iupred = np.empty((len(MSA), len(MSA[0])))

            for idx, record in enumerate(MSA):
                # Create map of ungapped to gapped coordinates
                idx_ungap = -1  # i == j at first non-gap
                bound_start_ungap = 0
                bound_start_gap = 0
                ungap_prev = record.seq[0] != '-'
                coords = []
                for idx_gap, sym in enumerate(record.seq):
                    ungap_curr = sym != '-'
                    idx_ungap += int(ungap_curr)
                    if ungap_curr is not ungap_prev:
                        gap_bound = (bound_start_gap, idx_gap)
                        ungap_bound = (bound_start_ungap, idx_ungap + 1) if ungap_prev else ()
                        coords.append((gap_bound, ungap_bound))

                        bound_start_ungap = idx_ungap
                        bound_start_gap = idx_gap
                    ungap_prev = ungap_curr
                gap_bound = (bound_start_gap, idx_gap + 1)
                ungap_bound = (bound_start_ungap, idx_ungap + 1) if ungap_prev else ()
                coords.append((gap_bound, ungap_bound))

                # Remove gaps from sequence
                seq = str(record.seq).translate({ord('-'): None})
                with tempfile.NamedTemporaryFile(mode='wt') as temp:
                    # Write sequence to tempfile
                    temp.write(seq)
                    temp.seek(0)  # Set stream to file start

                    # Execute IUPRED2a with tempfile
                    process = subprocess.run(['python', '../../../bin/iupred2a/iupred2a.py', temp.name, 'long'],
                                             capture_output=True, text=True, check=True)

                # Extract scores from output
                scores = []
                for line in process.stdout.rstrip().split('\n'):  # Remove trailing newline to prevent empty line
                    if line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    scores.append(float(fields[2]))
                scores = ndimage.gaussian_filter1d(scores, 2)

                # Fill IUPRED array with bools reflecting scores
                for coord_gap, coord_ungap in coords:
                    if coord_ungap == ():
                        MSA_iupred[idx, coord_gap[0]:coord_gap[1]] = -1
                    else:
                        MSA_iupred[idx, coord_gap[0]:coord_gap[1]] = scores[coord_ungap[0]:coord_ungap[1]]
            scores = np.any(np.logical_and(0 < MSA_iupred, MSA_iupred < 0.5), axis=0)

            # Find bounds
            bounds = []
            ordered_prev = scores[0]
            bound_start = 0
            for i, score in enumerate(scores):
                ordered_curr = score
                if ordered_curr is not ordered_prev:
                    bounds.append(((bound_start, i), ordered_prev))
                    bound_start = i
                ordered_prev = ordered_curr
            bounds.append(((bound_start, i + 1), ordered_curr))  # Bound for final segment

            # Create dataframe rows
            for bound, ordered in bounds:
                for record in MSA:
                    block_data.append({'ali_id': ali_id, 'seq_id': record.id, 'seg_id': hex(seg_num)[2:].zfill(8),
                                       'block_id': hex(block_num)[2:].zfill(8), 'bound': bound, 'ordered': ordered,
                                       'seq': record.seq[bound[0]:bound[1]]})
                    seg_num += 1
                block_num += 1

if not os.path.exists('out/'):
    os.mkdir('out/')

df = pd.DataFrame(block_data)
df.to_csv('out/segment_union.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../EggNOGv5_validation/filter_count/filter_count.py
    ../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv
../../EggNOGv5_validation/filter_unknown_realign/filter_unknown_realign.py
    ../../EggNOGv5_validation/filter_unknown_realign/out/align/
"""