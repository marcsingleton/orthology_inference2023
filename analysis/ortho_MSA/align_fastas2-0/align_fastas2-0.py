"""Align FASTAs of remaining OGs using longest isoforms."""

import multiprocessing as mp
import os
from subprocess import run, CalledProcessError
from time import time_ns


def run_cmd(file_id):
    cmd = (f'../../../bin/clustalo -i ../make_fastas2-0/out/{file_id}.tfa '
           f'--iterations 1 --threads=1 -o out/{file_id}.mfa')
    try:
        t0 = time_ns()
        run(cmd, shell=True, check=True)
        t1 = time_ns()
        return file_id, str(t1-t0)
    except CalledProcessError:
        return file_id, 'NaN'


num_processes = int(os.environ['SLURM_NTASKS'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        file_ids = [file[:-4] for file in os.listdir('../make_fastas2-0/out/')]
        rows = pool.map(run_cmd, file_ids)

    with open('out/times.tsv', 'w') as file:
        file.write('file_id\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../make_fastas2-0/make_fastas2-0.py
    ../make_fastas2-0/out/*.tfa
"""