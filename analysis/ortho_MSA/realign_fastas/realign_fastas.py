"""Align FASTAs of OGs."""

import multiprocessing as mp
import os
import re
from subprocess import run, CalledProcessError

import numpy as np
import scipy.ndimage as ndimage
import skbio
from src.utils import get_brownian_weights, read_fasta


def run_cmd(OGid):
    error_flag1 = False
    error_flag2 = False

    # Extract sequences and create initial alignment
    with open(f'out/{OGid}_temp1.fa', 'w') as file:
        for header, seq in read_fasta(f'../get_repseqs/out/{OGid}.afa'):
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            file.write(f'{header}\n{seqstring}\n')

    cmd = (f'../../../bin/mafft --globalpair --maxiterate 1000 '
           f'--thread 1 --anysymbol --allowshift --unalignlevel {unalignlevel1} --leavegappyregion '
           f'out/{OGid}_temp1.fa '
           f'1> out/{OGid}_temp1.afa 2> out/{OGid}_temp1.err')

    try:
        run(cmd, shell=True, check=True)
        input_msa = []
        for header, seq in read_fasta(f'out/{OGid}_temp1.afa'):
            ppid = re.search(ppid_regex, header).group(1)
            spid = re.search(spid_regex, header).group(1)
            input_msa.append({'header': header, 'ppid': ppid, 'spid': spid, 'seq': seq})
    except CalledProcessError:
        error_flag1 = True
        return OGid, error_flag1, error_flag2
    finally:
        os.remove(f'out/{OGid}_temp1.fa')
        os.remove(f'out/{OGid}_temp1.afa')
        os.remove(f'out/{OGid}_temp1.err')

    # Calculate weights
    spid_set = set([record['spid'] for record in input_msa])
    spid2ppids = {spid: [] for spid in spid_set}
    for record in input_msa:
        ppid, spid = record['ppid'], record['spid']
        spid2ppids[spid].append(ppid)

    tree = tree_template.shear(spid_set)
    spid_weights = {tip.name: weight for tip, weight in get_brownian_weights(tree)}
    ppid_weights = {}
    for spid, ppids in spid2ppids.items():
        weight = spid_weights[spid] / len(ppids)  # Species weight is distributed evenly over all associated proteins
        for ppid in ppids:
            ppid_weights[ppid] = weight

    # Create weight MSA
    weight_msa = np.zeros((len(input_msa), len(input_msa[0]['seq'])))
    for i, record in enumerate(input_msa):
        ppid, seq = record['ppid'], record['seq']
        weight = ppid_weights[ppid]
        for j, sym in enumerate(seq):
            if sym != '-':
                weight_msa[i, j] = weight

    # Extract poorly aligned regions
    regions = np.ones(len(input_msa[0]['seq']))
    binary = ndimage.binary_closing(weight_msa.sum(axis=0) >= threshold, structure=structure)  # Close on well-aligned regions
    for s, in ndimage.find_objects(ndimage.label(binary)[0]):
        if s.stop - s.start >= min_length:
            regions[s] = 0
    slices = [s for s, in ndimage.find_objects(ndimage.label(regions)[0])]

    # Align subsequences in regions and stitch together results
    idx = 0
    output_msa = {record['header']: [] for record in input_msa}
    for s in slices:
        # Append any aligned regions
        if idx < s.start:
            for record in input_msa:
                header, seq = record['header'], record['seq']
                output_msa[header].append(seq[idx:s.start])

        # Collect subsequences
        subseqs = []
        for record in input_msa:
            header, seq = record['header'], record['seq']
            subseq = seq[s]
            for sym in subseq:
                if sym != '-':
                    subseqs.append((header, subseq.replace('-', '')))
                    break

        # Align subsequences
        if len(subseqs) == 1:
            output = [(header, subseq) for header, subseq in subseqs]
        else:
            with open(f'out/{OGid}_temp2.fa', 'w') as file:
                for header, subseq in subseqs:
                    seqstring = '\n'.join([subseq[i:i+80] for i in range(0, len(subseq), 80)])
                    file.write(f'{header}\n{seqstring}\n')

            cmd = (f'../../../bin/mafft --globalpair --maxiterate 1000 '
                   f'--thread 1 --anysymbol --allowshift --unalignlevel {unalignlevel2} --leavegappyregion '
                   f'out/{OGid}_temp2.fa '
                   f'1> out/{OGid}_temp2.afa 2> out/{OGid}_temp2.err')
            try:
                run(cmd, shell=True, check=True)
                output = [(header, subseq) for header, subseq in read_fasta(f'out/{OGid}_temp2.afa')]
            except CalledProcessError:  # Catch re-alignment errors and use unaligned subseq; see notes for more details
                output = read_fasta(f'out/{OGid}_temp2.fa')
                error_flag2 = True
            finally:
                os.remove(f'out/{OGid}_temp2.fa')
                os.remove(f'out/{OGid}_temp2.afa')
                os.remove(f'out/{OGid}_temp2.err')

        length = len(output[0][1])
        output = {header: subseq for header, subseq in output}

        # Append aligned subsequences, filling missing subsequences with gaps
        for header, subseqs in output_msa.items():
            if header in output:
                subseq = output[header]
            else:
                subseq = length * '-'
            subseqs.append(subseq)

        idx = s.stop
    if idx < len(input_msa[0]['seq']):
        for record in input_msa:
            header, seq = record['header'], record['seq']
            output_msa[header].append(seq[idx:len(input_msa[0]['seq'])])

    # Write to file
    with open(f'out/{OGid}.afa', 'w') as file:
        for record in input_msa:  # Use same order as input
            header = record['header']
            seq = ''.join(output_msa[header])
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            file.write(f'{header}\n{seqstring}\n')

    return OGid, error_flag1, error_flag2


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

unalignlevel1 = 0.7
unalignlevel2 = 0.4
threshold = 0.5
structure = np.ones(3)
min_length = 10

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)

if __name__ == '__main__':
    OGids = []
    with open('../OG_filter/out/OG_filter.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            OGids.append(fields['OGid'])

    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        records = pool.map(run_cmd, OGids)

    with open('out/errors.tsv', 'w') as file:
        file.write('OGid\terror_flag1\terror_flag2\n')
        for record in records:
            file.write('\t'.join([str(field) for field in record]) + '\n')

"""
DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../make_fastas/make_fastas.py
    ../make_fastas/out/*.fa
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
"""