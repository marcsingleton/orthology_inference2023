""""Functions to calculate features associated with IDRs"""

import re
from math import log2

from localcider.sequenceParameters import SequenceParameters
from ipc import predict_isoelectric_point


# General functions
def count_group(seq, group):
    """Return count of residues matching amino acids in group."""
    count = 0
    for sym in seq:
        if sym in group:
            count += 1
    return count


def fraction_X(seq, x):
    """Return fraction of sequence matching amino acid x."""
    return seq.count(x) / len(seq)


def fraction_group(seq, group):
    """Return fraction of sequence matching amino acids in group."""
    count = count_group(seq, group)
    return count / len(seq)


def fraction_regex(seq, regex):
    """Return fraction of sequence matching a regular expression."""
    matches = re.findall(regex, seq)
    count = 0
    for match in matches:
        count += len(match)
    return count / len(seq)


# Amino acid content
def get_features_aa(seq):
    """Return fractions of sequence matching individual amino acids SPTAHQNG as dictionary."""
    features = {}
    for sym in 'SPTAHQNG':
        features['fraction_' + sym] = fraction_X(seq, sym)
    return features


# Charge properties
def count_positive(seq):
    """Return count of positively charged residues."""
    return count_group(seq, 'RK')


def count_negative(seq):
    """Return count of negatively charged residues."""
    return count_group(seq, 'DE')


def FCR(seq):
    """Return fraction of charged residues in sequence."""
    return (count_positive(seq) + count_negative(seq)) / len(seq)


def NCPR(seq):
    """Return net charge per residue in sequence."""
    return (count_positive(seq) - count_negative(seq)) / len(seq)


def net_charge(seq):
    """Return net charge of sequence."""
    return count_positive(seq) - count_negative(seq)


def net_charge_P(seq):
    """Return net charging accounting for likely phosphorylated serines and threonines."""
    psites = re.findall('[ST]P', seq)
    return net_charge(seq) - 1.5 * len(psites)


def RK_ratio(seq):
    """Return adjusted ratio of arginine to lysine residues."""
    r = 1 + seq.count('R')
    k = 1 + seq.count('K')
    return r / k


def ED_ratio(seq):
    """Return adjusted ratio of aspartate to glutamate residues."""
    e = 1 + seq.count('E')
    d = 1 + seq.count('D')
    return e / d


def get_features_charge(seq):
    """Return dictionary of all features associated with charge."""
    SeqOb = SequenceParameters(seq)
    return {'FCR': FCR(seq), 'NCPR': NCPR(seq),
            'net_charge': net_charge(seq), 'net_charge_P': net_charge_P(seq),
            'RK_ratio': RK_ratio(seq), 'ED_ratio': ED_ratio(seq),
            'kappa': SeqOb.get_kappa(), 'omega': SeqOb.get_Omega(), 'SCD': SeqOb.get_SCD()}


# Physiochemical properties
def fraction_acidic(seq):
    """Return fraction of acidic residues in sequence."""
    return fraction_group(seq, set('DE'))


def fraction_basic(seq):
    """Return fraction of basic residues in sequence."""
    return fraction_group(seq, set('RK'))


def fraction_aliphatic(seq):
    """Return fraction of aliphatic residues in sequence."""
    return fraction_group(seq, set('ALMIV'))


def fraction_aromatic(seq):
    """Return fraction of aromatic residues in sequence."""
    return fraction_group(seq, set('FYW'))


def fraction_polar(seq):
    """Return fraction of polar residues in sequence."""
    return fraction_group(seq, set('QNSTCH'))


def fraction_disorder(seq):
    """Return fraction of disorder-promoting residues in sequence."""
    return fraction_group(seq, set('TAGRDHQKSEP'))


def fraction_chainexp(seq):
    """Return fraction of chain-expanding residues in sequence."""
    return fraction_group(seq, set('EDRKP'))


def get_features_physchem(seq):
    """Return dictionary of all features associated with physiochemical properties."""
    SeqOb = SequenceParameters(seq)
    return {'fraction_acidic': fraction_acidic(seq), 'fraction_basic': fraction_basic(seq),
            'fraction_aliphatic': fraction_aliphatic(seq), 'fraction_aromatic': fraction_aromatic(seq),
            'fraction_polar': fraction_polar(seq), 'fraction_disorder': fraction_disorder(seq), 'fraction_chainexp': fraction_chainexp(seq),
            'hydropathy': SeqOb.get_uversky_hydropathy(), 'isopoint': predict_isoelectric_point(seq),
            'loglen': log2(len(seq)), 'PPII_propensity': SeqOb.get_PPII_propensity()}


# Sequence complexity
def get_features_complexity(seq):
    """Return dictionary of all features associated with sequence complexity."""
    repeats = ['Q', 'N', 'S', 'G', 'E', 'D', 'K', 'R', 'P',
               'QN', 'RG', 'FG', 'SG', 'SR', 'KAP', 'PTS']
    features = {}
    for repeat in repeats:
        features['repeat_' + repeat] = fraction_regex(seq, f'[{repeat}]' + '{2,}')
    features['wf_complexity'] = SequenceParameters(seq).get_linear_complexity(blobLen=len(seq))[1][0]  # Returns a 2xN matrix containing the complexity vector and the corresponding residue positions distributed equally along the sequence

    return features


# Summary
def get_features(seq):
    """Return dictionary of all features."""
    return {**get_features_aa(seq), **get_features_charge(seq), **get_features_physchem(seq), **get_features_complexity(seq)}
