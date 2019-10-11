"""Segment alignments by IUPRED2a and combine by union."""

import gzip
import numpy as np
import pandas as pd
import subprocess
import tempfile
from Bio import AlignIO
from scipy import ndimage

path = '../../EggNOGv5_validation/filter_count/7214_noX_members/10_10_members.tsv'
dir = '../../EggNOGv5_validation/filter_unknown_realign/align/'
thresh = 0.5

block_data = []  # Block data with raw sequences
block_num = 0  # Counter for numbering blocks

with open(path) as file:
    for line in file:
        fields = line.split('\t')
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
                    proc = subprocess.run(['python', '../../../bin/iupred2a/iupred2a.py', temp.name, 'long'],
                                          capture_output=True, text=True, check=True)

                # Extract scores from output
                scores = []
                for line in proc.stdout.rstrip().split('\n'):  # Remove trailing newline to prevent empty line
                    if line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    scores.append(float(fields[2]))

                # Fill IUPRED array with scores
                for coord_gap, coord_ungap in coords:
                    if coord_ungap is ():
                        MSA_iupred[idx, coord_gap[0]:coord_gap[1]] = -1
                    else:
                        MSA_iupred[idx, coord_gap[0]:coord_gap[1]] = scores[coord_ungap[0]:coord_ungap[1]]

            # Average IUPRED scores
            num_nongap = np.count_nonzero(MSA_iupred != -1, axis=0)
            MSA_avg = 1 / num_nongap * (np.sum(MSA_iupred, axis=0) + len(MSA) - num_nongap)
            MSA_avg = ndimage.gaussian_filter1d(MSA_avg, 2)  # Smooth final average

            # Find bounds
            bounds = []
            ordered_prev = MSA_avg[0] < thresh
            bound_start = 0
            for i, score in enumerate(MSA_avg):
                ordered_curr = score < thresh
                if ordered_curr is not ordered_prev:
                    bounds.append(((bound_start, i), ordered_prev))
                    bound_start = i
                ordered_prev = ordered_curr
            bounds.append(((bound_start, i + 1), ordered_curr))  # Bound for final segment

            # Create dataframe rows
            for bound, ordered in bounds:
                for record in MSA:
                    block_data.append({'ali_id': ali_id, 'seq_id': record.id, 'block_id': hex(block_num)[2:].zfill(8),
                                       'bound': bound, 'ordered': ordered, 'seq': record.seq[bound[0]:bound[1]]})
                block_num += 1

df = pd.DataFrame(block_data)
df.to_csv('segment_avg.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../EggNOGv5_validation/filter_count/filter_count.py
    ../../EggNOGv5_validation/filter_count/7214_noX_members/10_10_members.tsv
../../EggNOGv5_validation/filter_unknown_realign/filter_unknown_realign.py
    ../../EggNOGv5_validation/filter_unknown_realign/align/
"""