"""Extract HSPs from BLAST results"""

import multiprocessing as mp
import os
import re
from itertools import groupby, permutations


def parse_file(query_spid, subject_spid):
    output_hsps = []
    nulls = []
    with open(f'../blast_search/out/{query_spid}/{subject_spid}.blast') as file:
        query_ppid, input_hsps = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line == '# BLASTP 2.10.1+\n' and query_ppid is not None:  # Only add if previous search returned no hits
                    nulls.append({'qppid': query_ppid, 'qgnid': query_gnid, 'qspid': query_spid, 'sspid': subject_spid})
                elif line.startswith('# Query:'):
                    query_ppid = re.search(ppid_regex[genomes[query_spid]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record HSPs
            while line and not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                subject_ppid = re.search(ppid_regex[genomes[subject_spid]], fields[1]).group(1)
                subject_gnid = ppid2gnid[subject_ppid]
                values = [query_ppid, query_gnid, query_spid,
                          subject_ppid, subject_gnid, subject_spid,
                          *fields[2:],
                          False, False, False]
                input_hsps.append({column: f(value) for (column, f), value in zip(columns.items(), values)})
                line = file.readline()

            # Add best from HSP list (and catch cases where the last search returned no HSPs)
            if input_hsps:
                output_hsps.extend(filter_hsps(input_hsps))
            else:
                nulls.append({'qppid': query_ppid, 'qgnid': query_gnid, 'qspid': query_spid, 'sspid': subject_spid})

            if line.startswith('# BLAST processed'):
                break
            query_ppid, input_hsps = None, []  # Signals current search was successfully recorded

    # Write HSPs and nulls to file
    if not os.path.exists(f'out/hsps/{query_spid}/'):
        os.makedirs(f'out/hsps/{query_spid}/')
    if not os.path.exists(f'out/nulls/{query_spid}/'):
        os.makedirs(f'out/nulls/{query_spid}/')

    with open(f'out/hsps/{query_spid}/{subject_spid}.tsv', 'w') as file:
        file.write('\t'.join(columns) + '\n')
        for hsp in output_hsps:
            file.write('\t'.join([str(hsp[column]) for column in columns]) + '\n')
    with open(f'out/nulls/{query_spid}/{subject_spid}.tsv', 'w') as file:
        file.write('\t'.join(['qppid', 'qgnid', 'qspid', 'sspid']) + '\n')
        for null in nulls:
            file.write('\t'.join([null[column] for column in ['qppid', 'qgnid', 'qspid', 'sspid']]) + '\n')


def filter_hsps(input_hsps):
    # Sort and group HSPs by sppid and sgnid then extract max bitscore within groups
    input_hsps = sorted(input_hsps, key=lambda x: (x['sppid'], x['sgnid']))
    groups = {key: list(group) for key, group in groupby(input_hsps, lambda x: (x['sppid'], x['sgnid']))}
    bitscores = sorted([(ppid, gnid, max([hsp['bitscore'] for hsp in hsps])) for (ppid, gnid), hsps in groups.items()],
                       key=lambda x: x[2], reverse=True)

    # Get target bitscore cutoff
    target_cutoff = 0
    gnids = set()
    for ppid, gnid, bitscore in bitscores:
        if bitscore == bitscores[0][2]:  # Check if equal to max, i.e. the first in the sorted list
            gnids.add(gnid)  # Add gnid to allowable list if bitscore is equal to max
        if gnid not in gnids:
            target_cutoff = bitscore  # Choose target cutoff as the maximum bitscore associated with the next best gene
            break

    # Identify groups that pass target bitscore filter
    keys = []
    for ppid, gnid, bitscore in bitscores:
        if bitscore <= target_cutoff:  # Stop recording groups once bitscore is lower than target cutoff
            break
        keys.append((ppid, gnid))

    # Identify index, disjoint, and compatible HSPs
    output_hsps = []
    for key in keys:
        group = sorted(groups[key], key=lambda x: x['bitscore'], reverse=True)
        group[0]['index_hsp'] = True  # Mark "best" HSP as index for subsequent filtering

        disjoint_hsps = []
        for hsp in group:
            if is_disjoint(hsp, disjoint_hsps, 'query') and is_disjoint(hsp, disjoint_hsps, 'subject'):
                hsp['disjoint'] = True  # Mark HSPs as disjoint greedily, beginning with highest score
                if hsp['bitscore'] >= compatible_cutoff:
                    hsp['compatible'] = True  # Disjoint HSPs are compatible if they pass the cutoff
                disjoint_hsps.append(hsp)

        compatible_hsps = [hsp for hsp in group if hsp['compatible']]
        for hsp in group:
            if is_compatible(hsp, compatible_hsps, 'query') and is_compatible(hsp, compatible_hsps, 'subject') and hsp['bitscore'] > compatible_cutoff:
                hsp['compatible'] = True
                compatible_hsps.append(hsp)

        output_hsps.extend(group)  # Add all HSPs in this group to output for query sequence

    return output_hsps


def is_disjoint(hsp, hsp_list, key_type):
    """Return if hsp is disjoint with all HSPs in hsp_list.

    For each test_hsp in hsp_list, hsp must start after test_hsp ends or end
    before test_hsp starts. If so, returns True, else returns False.
    """
    if key_type == 'query':
        start_key = 'qstart'
        end_key = 'qend'
    elif key_type == 'subject':
        start_key = 'sstart'
        end_key = 'send'
    else:
        raise ValueError('key_type is not query or subject')

    for test_hsp in hsp_list:
        if not (hsp[start_key] > test_hsp[end_key] or test_hsp[start_key] > hsp[end_key]):
            return False
    return True


def is_compatible(hsp, hsp_list, key_type):
    """Return if hsp is compatible with all HSPs in hsp_list.

    For each test_hsp in hsp_list, overlap of hsp and test_hsp must not exceed
    50% of the length of either. If so, returns True, else returns False.
    """
    if key_type == 'query':
        start_key = 'qstart'
        end_key = 'qend'
    elif key_type == 'subject':
        start_key = 'sstart'
        end_key = 'send'
    else:
        raise ValueError('key_type is not query or subject')

    for test_hsp in hsp_list:
        if not (hsp[start_key] > test_hsp[end_key] or test_hsp[start_key] > hsp[end_key]):
            start = max(hsp[start_key], test_hsp[start_key])
            end = min(hsp[end_key], test_hsp[end_key])
            length = end-start
            if length / (hsp[end_key]-hsp[start_key]) >= 0.5 or length / (test_hsp[end_key]-test_hsp[start_key]) >= 0.5:
                return False
    return True


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}
columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'length': int, 'nident': int, 'gaps': int,
           'qlen': int, 'qstart': int, 'qend': int,
           'slen': int, 'sstart': int, 'send': int,
           'evalue': float, 'bitscore': float,
           'index_hsp': bool, 'disjoint': bool, 'compatible': bool}
compatible_cutoff = 50  # Bitscore cutoff for accepting HSPs as compatible
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

# Load genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, _ = line.rstrip('\n').split('\t')
        genomes[spid] = source

# Load sequence data
ppid2gnid = {}
with open('../sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.rstrip('\n').split('\t')
        ppid2gnid[ppid] = gnid

# Parse BLAST results
if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(parse_file, permutations(genomes, 2))

"""
NOTES
This script makes a few decisions about what exactly a best hit to both a polypeptide and a gene are. Regarding
polypeptides, since it is possible for there to be multiple hits to different segments of the sequence, the "best" is
taken as the HSP with the largest bitscore. However, all other HSPs from this hit are recorded. Proceeding from highest
to lowest bitscore, these HSPs are marked as disjoint if they do not overlap with any HSPs previously marked as
disjoint. Additionally other polypeptide hits are recorded if the strongest HSP in that hit exceeds the bitscore of
the strongest HSP to a polypeptide that is not associated with the same gene as the best polypeptide hit. In other
words, this approach takes a broad "gene-centric" approach to initially filtering the raw BLAST output where all hits
associated with the best gene are recorded as long as the bitscore of their strongest HSP exceeds that of the next best
gene. If there are multiple HSPs with the same best bitscore, then they are both treated as the best hit. This means a
single polypeptide can have hits to multiple genes which in turn contain multiple polypeptides which in turn contain
multiple HSPs.

Regarding genes, a strict adherence to the top hit criterion would allow only polypeptide per gene. This can create
mutually exclusive sets of polypeptide hits within a single orthologous group, a more flexible strategy is allowing
multiple best hits within a gene. The most natural extension of the polypeptide criterion is choosing top hits ranked by
evalue until the parent gene of the polypeptide changes. Subsequent hits to the top gene are considered ambiguous and
ignored (similar to how subsequent hits to a polypeptide are ignored in the polypeptide case). This criterion will
exclude some polypeptides that are associated with the best gene, but since these are ranked lower than hits to
polypeptides in other genes, they should not be included since other stronger hits are discarded.

DEPENDENCIES
../config/genomes.tsv
../blast_search/blast_search.py
    ../blast_search/out/*
../sequence_data/sequence_data.py
    ../sequence_data/out/sequence_data.tsv
"""