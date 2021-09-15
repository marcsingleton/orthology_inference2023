"""Run AUCpreD on individual sequences in trimmed alignments."""

import multiprocessing as mp
import os
import re
import subprocess


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def run_aucpred(OGid):
    msa = load_msa(f'../insertion_trim/out/{OGid}.mfa')
    prefix = f'out/raw/{OGid}/'

    if not os.path.exists(prefix):
        os.mkdir(prefix)

    for header, seq in msa:
        ppid = re.search(r'ppid=([A-Za-z0-9_]+)', header).group(1)
        seq = seq.translate({ord('-'): None, ord('.'): None})
        if len(seq) < 10000:  # AUCpreD uses PSIPRED which has a length limit of 10000
            with open(f'{prefix}{ppid}.fasta', 'w') as file:
                seqstring = '\n'.join([seq[i:i + 80] for i in range(0, len(seq), 80)]) + '\n'
                file.write(header + '\n' + seqstring)
            subprocess.run(f'../../../bin/Predict_Property/AUCpreD.sh -i {prefix}{ppid}.fasta -o {prefix}',
                           check=True, shell=True)
            os.remove(f'{prefix}{ppid}.fasta')


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if __name__ == '__main__':
    if not os.path.exists('out/raw/'):
        os.makedirs('out/raw/')

    with mp.Pool(processes=num_processes) as pool:
        OGids = [path.split('.')[0] for path in os.listdir('../insertion_trim/out/') if path.endswith('.mfa')]
        pool.map(run_aucpred, OGids)

"""
DEPENDENCIES
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""