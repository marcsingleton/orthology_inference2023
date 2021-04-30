"""Segment alignments by IUPRED2a score and write subsequences to dataframe."""

import gzip
import os
import pandas as pd
import subprocess
import tempfile
from Bio import AlignIO
from scipy import ndimage

path = '../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv'
dir = '../../EggNOGv5_validation/filter_unknown_realign/out/align/'
thresh = 0.5

seg_data = []  # Segment data with raw sequences
seg_num = 0  # Counter for numbering rows

with open(path) as file:
    for line in file:
        fields = line.split('\t')
        ali_id = fields[1]
        with gzip.open(dir + ali_id + '.raw_alg.faa.gz', 'rt') as file_MSA:  # read, text mode
            MSA = AlignIO.read(file_MSA, 'fasta')
            for record in MSA:
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

                # Find bounds
                bounds = []
                ordered_prev = scores[0] < thresh
                bound_start = 0
                for i, score in enumerate(scores):
                    ordered_curr = score < thresh
                    if ordered_curr is not ordered_prev:
                        bounds.append(((bound_start, i), ordered_prev))
                        bound_start = i
                    ordered_prev = ordered_curr
                bounds.append(((bound_start, i + 1), ordered_curr))  # Bound for final segment

                # Create dataframe rows
                for bound, ordered in bounds:
                    seg_data.append({'ali_id': ali_id, 'seq_id': record.id, 'seg_id': hex(seg_num)[2:].zfill(8),
                                     'bound': bound, 'ordered': ordered, 'seq': seq[bound[0]:bound[1]]})
                    seg_num += 1

if not os.path.exists('out/'):
    os.mkdir('out/')

df = pd.DataFrame(seg_data)
df.to_csv('segment_iupred2a.tsv', sep='\t', index=False)

"""
NOTES
Must run in bash shell to allow correct behavior of NamedTemporaryFile

DEPENDENCIES
../../EggNOGv5_validation/filter_count/filter_count.py
    ../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv
../../EggNOGv5_validation/filter_unknown_realign/filter_unknown_realign.py
    ../../EggNOGv5_validation/filter_unknown_realign/out/align/
"""