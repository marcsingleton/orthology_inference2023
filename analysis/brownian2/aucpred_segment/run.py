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
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def run_aucpred(OGid):
    msa = load_msa(f'../trim_extract/out/{OGid}.mfa')

    if not os.path.exists(f'out/{OGid}/'):
        os.mkdir(f'out/{OGid}/')

    for header, seq in msa:
        ppid = re.match(r'>ppid=([A-Za-z0-9_]+)\|', header).group(1)
        seq = seq.translate({ord('-'): None, ord('.'): None})
        if len(seq) < 10000:  # AUCpreD uses PSIPRED which has a length limit of 10000
            with open(f'out/{OGid}/{ppid}.fasta', 'w') as file:
                seqstring = '\n'.join([seq[i:i + 80] for i in range(0, len(seq), 80)]) + '\n'
                file.write(header + seqstring)
            subprocess.run(f'../../../bin/Predict_Property/AUCpreD.sh -i out/{OGid}/{ppid}.fasta -o out/{OGid}/',
                           check=True, shell=True)
            os.remove(f'out/{OGid}/{ppid}.fasta')


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        OGids = [path.split('.')[0] for path in os.listdir('../trim_extract/out/') if path.endswith('.mfa')]
        pool.map(run_aucpred, OGids)

"""
DEPENDENCIES
../trim_extract/extract.py
    ../trim_extract/out/*.mfa
"""