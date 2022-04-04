"""Extract HSPs from BLAST results"""

import multiprocessing as mp
import os
import re
from itertools import groupby, permutations


def parse_file(query_spid, subject_spid):
    output_hsps, nulls = [], []
    with open(f'../blast_search/out/{query_spid}/{subject_spid}.blast') as file:
        query_ppid, input_hsps = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                is_null = (line == f'# BLASTP {blast_version}\n' or line.startswith('# BLAST processed')) and query_ppid is not None
                if is_null:  # Only True for two consecutive headers or header followed by file end (# BLAST processed)
                    nulls.append({'qppid': query_ppid, 'qgnid': query_gnid})
                elif line.startswith('# Query:'):
                    query_ppid = re.search(ppid_regex[genomes[query_spid]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record HSPs
            while line and not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                subject_ppid = ppid2ppid[fields[1]]
                subject_gnid = ppid2gnid[subject_ppid]
                values = [query_ppid, query_gnid,
                          subject_ppid, subject_gnid,
                          *fields[2:],
                          False, False, False]
                input_hsps.append({column: f(value) for (column, f), value in zip(columns.items(), values)})
                line = file.readline()
            output_hsps.extend(filter_hsps(input_hsps))

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
        file.write('\t'.join(['qppid', 'qgnid']) + '\n')
        for null in nulls:
            file.write('\t'.join([null[column] for column in ['qppid', 'qgnid']]) + '\n')


def filter_hsps(input_hsps):
    # Sort and group HSPs by sppid and sgnid then extract max bitscore within groups
    hsp2key = lambda x: (x['sppid'], x['sgnid'])
    input_hsps = sorted(input_hsps, key=hsp2key)
    grouped_hsps = {key: list(group) for key, group in groupby(input_hsps, hsp2key)}

    # Create records for identifying bitscore cutoff
    records = []
    for (ppid, gnid), hsps in grouped_hsps.items():
        bitscore = max([hsp['bitscore'] for hsp in hsps])
        records.append((ppid, gnid, bitscore))
    records = sorted(records, key=lambda x: x[2], reverse=True)
    max_bitscore = records[0][2]

    # Get subject bitscore cutoff
    # (This is a separate step to account for cases where a PPID associated with an allowable GNID has a bitscore equal to the cutoff)
    subject_cutoff, gnids = 0, set()
    for ppid, gnid, bitscore in records:
        if bitscore == max_bitscore:
            gnids.add(gnid)  # Add gnid to allowable list if bitscore is equal to max
        if gnid not in gnids:
            subject_cutoff = bitscore  # Choose subject cutoff as the maximum bitscore associated with the next best gene
            break

    # Identify groups that pass subject bitscore filter
    keys = []
    for ppid, gnid, bitscore in records:
        if bitscore <= subject_cutoff:  # Stop recording groups once bitscore is lower than subject cutoff
            break
        keys.append((ppid, gnid))

    # Identify index, disjoint, and compatible HSPs
    output_hsps = []
    for key in keys:
        group = sorted(grouped_hsps[key], key=lambda x: x['bitscore'], reverse=True)
        group[0]['index_hsp'] = True  # Mark "best" HSP as index for subsequent filtering

        disjoint_hsps = []
        for hsp in group:
            if is_disjoint(hsp, disjoint_hsps, 'query') and is_disjoint(hsp, disjoint_hsps, 'subject'):
                hsp['disjoint'] = True  # Mark HSPs as disjoint greedily, beginning with highest score
                if hsp['evalue'] <= compatible_cutoff:
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
              'NCBI': r'>([NXY]P_[0-9]+(\.[0-9]+)?)'}
columns = {'qppid': str, 'qgnid': str,
           'sppid': str, 'sgnid': str,
           'length': int, 'nident': int, 'gaps': int,
           'qlen': int, 'qstart': int, 'qend': int,
           'slen': int, 'sstart': int, 'send': int,
           'evalue': float, 'bitscore': float,
           'index_hsp': bool, 'disjoint': bool, 'compatible': bool}
blast_version = '2.13.0+'  # Used to identify headers of queries
compatible_cutoff = 1E-10  # E-value cutoff for accepting HSPs as compatible
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

# Load genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, _ = line.rstrip('\n').split('\t')
        genomes[spid] = source

# Load sequence data
ppid2gnid, ppid2ppid = {}, {}
with open('../sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.rstrip('\n').split('\t')
        ppid2gnid[ppid] = gnid
        ppid2ppid[ppid.split('.')[0]] = ppid  # Mapping from PPID w/o version to PPID w/ version since truncated in BLAST output

# Parse BLAST results
if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(parse_file, permutations(genomes, 2))

"""
NOTES
The main idea for identifying orthologs is to use clusters of reciprocal best hits from pairwise database searches of
the genomes. Because BLAST returns local alignments (i.e. alignments between portions of the query and subject sequences
rather than alignments from end-to-end), the best hit for a query sequence is not immediately obvious since it is
potentially split across multiple local alignments, called HSPs or high-scoring segment pairs. Additionally, since hits
must overlap by 50% to be included, not accounting for multiple HSPs could miss some highly significant hits.

To reduce the workload of merging HSPs into hits later on, the raw BLAST results are first filtered into a reduced list
of HSPs which are likely to contain the best hit. The strategy is gene- and HSP-centric, where the maximum bitscore of
the HSP associated with the next best gene sets the cutoff. As long as a group of HSPs has one which exceeds that
cutoff, they are all included for the next step of processing. This means that no single isoform associated with the
best gene can set the cutoff, but some isoforms may not be included if their best HSP is weaker

DEPENDENCIES
../config/genomes.tsv
../blast_search/blast_search.py
    ../blast_search/out/*
../sequence_data/sequence_data.py
    ../sequence_data/out/sequence_data.tsv
"""