""""Features set used in feature_calc.py."""

import re
from localcider.sequenceParameters import SequenceParameters
from math import log2
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# General functions
def remove_gaps(seq):
    return seq.translate({ord('-'): None})


def count_group(seq, group):
    count = 0
    for sym in group:
        count += seq.count(sym)
    return count


def frac_X(seq, x):
    """Return fraction of sequence for arbitrary amino acid X."""
    return seq.count(x) / len(seq)


def frac_group(seq, group):
    """Return fraction of sequence for arbitrary group of amino acids where the group is a string of amino acid symbols."""
    count = 0
    for sym in group:
        count += seq.count(sym)
    return count / len(seq)


def frac_pattern(seq, pattern):
    """Return fraction of sequence matching a regular expression."""
    matches = re.findall(pattern, seq)
    sum = 0
    for match in matches:
     sum += len(match)
    return sum / len(seq)


# Amino acid content
def frac_aa(seq):
    """Return fractions of SPTAHQNG as dictionary."""
    fracs = {}
    for sym in 'SPTAHQNG':
        fracs['frac_' + sym] = frac_X(seq, sym)
    return fracs


# Charge properties
def count_pos(seq):
    return count_group(seq, 'RK')


def count_neg(seq):
    return count_group(seq, 'DE')


def FCR(seq):
    return (count_pos(seq) + count_neg(seq)) / len(seq)


def NCPR(seq):
    return (count_pos(seq) - count_neg(seq)) / len(seq)


def net_charge(seq):
    return count_pos(seq) - count_neg(seq)


def net_charge_P(seq):
    psites = re.findall('[ST]P', seq)
    return net_charge(seq) - 1.5 * len(psites)


def RK_ratio(seq):
    r = 1 + seq.count('R')
    k = 1 + seq.count('K')
    return r / k


def ED_ratio(seq):
    e = 1 + seq.count('E')
    d = 1 + seq.count('D')
    return e / d


def feat_charge(seq):
    SeqOb = SequenceParameters(seq)
    return {'FCR': FCR(seq), 'NCPR': NCPR(seq), 'net_charge': net_charge(seq), 'net_charge_P': net_charge_P(seq),
            'RK_ratio': RK_ratio(seq), 'ED_ratio': ED_ratio(seq), 'kappa': SeqOb.get_kappa(), 'omega': SeqOb.get_Omega(),
            'SCD': SeqOb.get_SCD()}


# Physiochemical properties
def frac_acidic(seq):
    return frac_group(seq, 'DE')


def frac_basic(seq):
    return frac_group(seq, 'RK')


def frac_aliphatic(seq):
    return frac_group(seq, 'ALMIV')


def frac_chainexp(seq):
    return frac_group(seq, 'EDRKP')


def frac_polar(seq):
    return frac_group(seq, 'QNSTCH')


def frac_aromatic(seq):
    return frac_group(seq, 'FYW')


def frac_disorder(seq):
    return frac_group(seq, 'TAGRDHQKSEP')


def feat_physchem(seq):
    SeqOb = SequenceParameters(seq)
    return {'frac_acidic': frac_acidic(seq), 'frac_basic': frac_basic(seq),
            'frac_aliphatic': frac_aliphatic(seq), 'frac_chainexp': frac_chainexp(seq),
            'frac_polar': frac_polar(seq), 'frac_aromatic': frac_aromatic(seq),
            'frac_disorder': frac_disorder(seq), 'loglen': log2(len(seq)),
            'hydropathy': SeqOb.get_uversky_hydropathy(), 'iso_point': ProteinAnalysis(seq).isoelectric_point(),
            'PPII_prop': SeqOb.get_PPII_propensity()}


# Repeats and complexity
def rep_fractions(seq):
    """Return fractions of repeats as dictionary."""
    fracs = {}
    repeats = ['Q', 'N', 'S', 'G', 'E', 'D', 'K', 'R', 'P',
               'QN', 'RG', 'FG', 'SG', 'SR', 'KAP', 'PTS']
    for repeat in repeats:
        fracs['rep_' + repeat] = frac_pattern(seq, f'[{repeat}]' + '{2,}')
    return fracs


def feat_complexity(seq):
    return {'wf_complexity': SequenceParameters(seq).get_linear_complexity(blobLen=len(seq))[1][0],  # Returns a 2xN matrix containing the complexity vector and the corresponding residue positions distributed equally along the sequence
            **rep_fractions(seq)}


# Summary
def feat_all(seq):
    return {**frac_aa(seq), **feat_charge(seq), **feat_physchem(seq), **feat_complexity(seq)}
