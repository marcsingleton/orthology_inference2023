"""Add in-paralogs to clusters."""

import os
import re
from itertools import groupby, product

import numpy as np


def parse_file(query_spid, subject_spid):
    hits = []
    with open(f'../../ortho_search/blast_search/out/{query_spid}/{subject_spid}.blast') as file:
        query_ppid, input_hsps = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line.startswith('# Query:'):
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
            if input_hsps:  # Guard against empty input which can occur on final pass through block at file end
                hits.extend(filter_hsps(input_hsps))

            query_ppid, input_hsps = None, []  # Signals current search was successfully recorded

    return hits


def filter_hsps(input_hsps):
    # Sort and group HSPs by sppid and sgnid then extract max bitscore within groups
    hsp2key = lambda x: (x['sppid'], x['sgnid'])
    input_hsps = sorted([hsp for hsp in input_hsps if hsp['qgnid'] != hsp['sgnid']], key=hsp2key)  # Remove self HSPs
    grouped_hsps = {key: list(group) for key, group in groupby(input_hsps, hsp2key)}

    # Identify index, disjoint, and compatible HSPs
    hits = []
    for group in grouped_hsps.values():
        group = sorted(group, key=lambda x: x['bitscore'], reverse=True)
        group[0]['index_hsp'] = True  # Mark "best" HSP as index for subsequent filtering

        disjoint_hsps = []
        for hsp in group:
            if is_disjoint(hsp, disjoint_hsps, 'query') and is_disjoint(hsp, disjoint_hsps, 'subject') and hsp['evalue'] <= evalue_cutoff:
                hsp['disjoint'] = True  # Mark HSPs as disjoint greedily, beginning with highest score
                hsp['compatible'] = True  # All disjoint HSPs are compatible
                disjoint_hsps.append(hsp)

        compatible_hsps = [hsp for hsp in group if hsp['compatible']]
        for hsp in group:
            if is_compatible(hsp, compatible_hsps, 'query') and is_compatible(hsp, compatible_hsps, 'subject') and hsp['evalue'] <= evalue_cutoff:
                hsp['compatible'] = True
                compatible_hsps.append(hsp)

        if compatible_hsps:  # Check for empty list
            hits.append(hsps2hit(compatible_hsps))

    return hits


def hsps2hit(hsps):
    # Calculate values from all HSPs
    hit = {key: hsps[0][key] for key in ['qppid', 'qgnid', 'sppid', 'sgnid', 'qlen', 'slen']}
    hit['chspnum'] = len(hsps)
    qcov = np.zeros((1, hsps[0]['qlen']), dtype=bool)
    scov = np.zeros((1, hsps[0]['slen']), dtype=bool)
    for hsp in hsps:
        qcov[0, hsp['qstart']-1:hsp['qend']] = True
        scov[0, hsp['sstart']-1:hsp['send']] = True
    hit['cnqa'] = qcov.sum()
    hit['cnsa'] = scov.sum()

    # Calculate values from disjoint HSPs only
    disjoint_hsps = [hsp for hsp in hsps if hsp['disjoint']]
    hit['hspnum'] = len(disjoint_hsps)
    hit['nqa'] = sum([hsp['qend'] - hsp['qstart'] + 1 for hsp in disjoint_hsps])
    hit['nsa'] = sum([hsp['send'] - hsp['sstart'] + 1 for hsp in disjoint_hsps])
    hit['bitscore'] = sum([hsp['bitscore'] for hsp in disjoint_hsps])

    return hit


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


def is_reciprocal(qppid, sppid, graph):
    # Check both directions since not all hits are in graph
    try:
        reciprocal1 = qppid in graph[sppid]
        reciprocal2 = sppid in graph[qppid]
        reciprocal = reciprocal1 and reciprocal2
    except KeyError:
        reciprocal = False
    return reciprocal


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}
columns = {'qppid': str, 'qgnid': str,
           'sppid': str, 'sgnid': str,
           'length': int, 'nident': int, 'gaps': int,
           'qlen': int, 'qstart': int, 'qend': int,
           'slen': int, 'sstart': int, 'send': int,
           'evalue': float, 'bitscore': float,
           'index_hsp': bool, 'disjoint': bool, 'compatible': bool}
blast_version = '2.13.0+'  # Used to identify headers of queries
evalue_cutoff = 1E-10  # E-value cutoff for accepting HSPs as disjoint or compatible

# Load genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, _ = line.rstrip('\n').split('\t')
        genomes[spid] = source

# Load sequence data
ppid2gnid, ppid2ppid = {}, {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, _ = line.rstrip('\n').split('\t')
        ppid2gnid[ppid] = gnid
        ppid2ppid[ppid.split('.')[0]] = ppid  # Mapping from PPID w/o version to PPID w/ version since truncated in BLAST output

# Load best scores from graph
ppid2score = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ppid2score[node] = max([float(adj.split(':')[1]) for adj in adjs.split(',')])

# Load OGs
OGs = []
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, algorithm, edges = line.rstrip('\n').split('\t')
        edges = [edge.split(':') for edge in edges.split(',')]
        OGs.append((component_id, OGid, algorithm, edges))

# Parse raw BLAST output finding intra-genome hits with scores greater than best inter-genome hit
graph = {}
for spid in genomes:
    hits = parse_file(spid, spid)
    for (qppid, _), group in groupby(hits, lambda x: (x['qppid'], x['sgnid'])):  # Group by gene as well to select best isoforms only
        group = list(group)  # Convert from iterator to list since used multiple times
        max_bitscore = max([hit['bitscore'] for hit in group])
        if max_bitscore <= ppid2score.get(qppid, 0):  # Only add in-paralogs to graph if they exceed score of inter-genome hits
            continue
        for hit in group:
            if max_bitscore == hit['bitscore']:  # Only record hits with maximum bitscore
                qppid, sppid = hit['qppid'], hit['sppid']
                qlen, cnqa = hit['qlen'], hit['cnqa']
                bitscore = hit['bitscore']

                if cnqa / qlen >= 0.5:
                    try:
                        graph[qppid].add(sppid)
                    except KeyError:
                        graph[qppid] = {sppid}

# Find in-paralogs
ppid2paralogs = {}
for node, adjs in graph.items():
    ppid2paralogs[node] = [node]
    for adj in adjs:
        if is_reciprocal(node, adj, graph):
            ppid2paralogs[node].append(adj)

# Add in-paralogs to clusters
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/clusters.tsv', 'w') as file:
    file.write('component_id\tOGid\talgorithm\tedges\n')
    for component_id, OGid, algorithm, edges1 in OGs:
        edges2 = []
        for node1, node2 in edges1:
            paralogs1 = ppid2paralogs.get(node1, [node1])  # Not all PPIDs are in graph so return [self]
            paralogs2 = ppid2paralogs.get(node2, [node2])
            edges2.extend(product(paralogs1, paralogs2))
        edgestring = ','.join([f'{node1}:{node2}' for node1, node2 in edges2])
        file.write(f'{component_id}\t{OGid}\t{algorithm}\t{edgestring}\n')

"""
DEPENDENCIES
../../ortho_search/blast_search/blast_search.py
    ../../ortho_search/blast_search/out/*
../config/genomes.tsv
../cluster4+_graph/cluster.py
    ../cluster4+_graph/out/4clique/clusters.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
../sequence_data/sequence_data.py
    ../sequence_data/out/sequence_data.tsv
"""