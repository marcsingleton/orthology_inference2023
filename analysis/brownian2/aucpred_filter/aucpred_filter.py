"""Remove sequences and regions from segments that do not pass quality filters."""

import re
import os

import scipy.stats as stats


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


def spid_filter(spids):
    conditions = [({'dnov', 'dvir'}, 1),
                  ({'dmoj', 'dnav'}, 1),
                  ({'dinn', 'dgri', 'dhyd'}, 2),
                  ({'dgua', 'dsob'}, 1),
                  ({'dbip', 'dana'}, 1),
                  ({'dser', 'dkik'}, 1),
                  ({'dele', 'dfik'}, 1),
                  ({'dtak', 'dbia'}, 1),
                  ({'dsuz', 'dspu'}, 1),
                  ({'dsan', 'dyak'}, 1),
                  ({'dmel'}, 1),
                  ({'dmau', 'dsim', 'dsec'}, 1)]
    return all([len(spids & group) >= num for group, num in conditions])


min_length = 30
alpha = 0.001
ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'

# Load regions
OGid2regions = {}
with open('../aucpred_segment/out/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.split()
        try:
            OGid2regions[OGid].append((start, stop, disorder))
        except KeyError:
            OGid2regions[OGid] = [(start, stop, disorder)]

# Filter regions
records = []
for OGid, regions in OGid2regions.items():
    msa = load_msa(f'../insertion_trim/out/{OGid}.mfa')

    for region in regions:
        # Get indices and length
        start, stop = int(region[0]), int(region[1])
        length = stop - start
        disorder = region[2]

        # Extract segments and counts of gaps and (non-gap) symbols
        segments0 = []
        for header, seq in msa:
            ppid = re.search(ppid_regex, header).group(1)
            spid = re.search(spid_regex, header).group(1)
            segment = seq[start:stop]
            gaps = sum([segment.count(sym) for sym in ['-', '.']])
            syms = length - gaps

            segments0.append((ppid, spid, segment, gaps, syms))

        # 1 Filter by length
        segments1 = []
        for record in segments0:
            syms = record[4]
            if syms > min_length:
                segments1.append(record)
        if not segments1:
            continue

        # 2 Filter by gap composition (Pearson's chi-squared test)
        total_gaps = sum([record[3] for record in segments1]) + 1  # Pseudocounts prevent later division by zero
        total_syms = sum([record[4] for record in segments1]) + 1
        total = total_gaps + total_syms
        p_gaps = total_gaps / total
        p_syms = total_syms / total
        exp_gaps = p_gaps * length
        exp_syms = p_syms * length

        segments2 = []
        for ppid, spid, segment, gaps, syms in segments1:
            chi2 = (gaps - exp_gaps) ** 2 / exp_gaps + (syms - exp_syms) ** 2 / exp_syms
            pvalue = 1 - stats.chi2.cdf(chi2, df=1)
            if pvalue >= alpha:
                segments2.append((ppid, spid, segment))
            else:
                print(f'Rejected {ppid} in {OGid}:{start}:{stop} (p = {pvalue})!')
        if not segments2:
            continue

        # 3 Filter by phylogenetic diversity
        ppids = [ppid for ppid, _, _ in segments2]
        spids = set([spid for _, spid, _ in segments2])
        if len(spids) == 20 and spid_filter(spids):
            records.append((OGid, str(start), str(stop), str(disorder), ','.join(ppids)))

# Write records to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/segments.tsv', 'w') as file:
    fields = ['OGid', 'start', 'stop', 'disorder', 'ppids']
    file.write('\t'.join(fields) + '\n')
    for record in records:
        file.write('\t'.join(record) + '\n')

"""
DEPENDENCIES
../aucpred_segment/segment.py
    ../aucpred_segment/out/segments.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""