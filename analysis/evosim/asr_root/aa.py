"""Calculate ancestral sequence reconstruction at root for amino acid process."""

import os
import re
from io import StringIO

import numpy as np
import skbio
from scipy.special import gammainc
from scipy.stats import gamma
from src.evosim.asr import get_conditional, get_tree
from src.utils import read_fasta


def load_model(path):
    with open(path) as file:
        # Parse exchangeability matrix
        matrix = np.zeros((len(alphabet), len(alphabet)))
        for i in range(len(alphabet)-1):
            line = file.readline()
            for j, value in enumerate(line.split()):
                matrix[i + 1, j] = float(value)
                matrix[j, i + 1] = float(value)

        # Parse equilibrium frequencies
        for _ in range(2):
            line = file.readline()
        freqs = np.array([float(value) for value in line.split()])
    matrix = freqs * matrix
    rate = (freqs * matrix.sum(axis=1)).sum()
    matrix = matrix / rate  # Normalize average rate to 1
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return matrix, freqs


alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
sym2idx = {sym: idx for idx, sym in enumerate(alphabet)}

models = {'WAG': load_model('../config/WAG.txt'), '../config/50red_D.txt': load_model('../config/50red_D.txt')}

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path[:-7] for path in os.listdir('../asr_aa/out/') if path.endswith('.iqtree')]
for OGid in OGids:
    # Load trees
    tree1 = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
    with open(f'../asr_aa/out/{OGid}.iqtree') as file:
        line = file.readline()
        while line != 'Tree in newick format:\n':
            line = file.readline()
        for _ in range(2):
            line = file.readline()
    tree2 = skbio.read(StringIO(line), 'newick', skbio.TreeNode)
    tree = get_tree(tree1, tree2)

    # Load partitions
    partitions = {}

    # Load partition model parameters
    with open(f'../asr_aa/out/{OGid}.iqtree') as file:
        # Get partition ID and name
        line = file.readline()
        while line.split() != ['ID', 'Name', 'Type', 'Seq', 'Site', 'Unique', 'Infor', 'Invar', 'Const']:  # Spacing can differ so check for fields
            line = file.readline()
        line = file.readline()
        while line != '\n':
            fields = line.split()  # File is not explicitly delimited, so just split on whitespace
            partition_id, name = int(fields[0]), fields[1]
            partitions[partition_id] = {'name': name}
            line = file.readline()

        # Get partition model parameters
        while line != '  ID  Model           Speed  Parameters\n':
            line = file.readline()
        line = file.readline()
        while line != '\n':
            fields = line.split()  # File is not explicitly delimited, so just split on whitespace
            partition_id, speed, parameters = int(fields[0]), float(fields[2]), fields[3]
            match = re.search(r'(?P<model>[^+]+)\+I{(?P<pinv>[0-9.e-]+)}\+G(?P<num_categories>[0-9]+){(?P<alpha>[0-9.e-]+)}', parameters)
            partition = partitions[partition_id]
            partition.update({'model': match['model'], 'speed': speed,
                              'pinv': float(match['pinv']), 'alpha': float(match['alpha']), 'num_categories': int(match['num_categories'])})
            line = file.readline()

    # Load partition regions
    with open(f'../asr_aa/out/{OGid}.nex') as file:
        partition_id = 1
        for line in file:
            if 'charset' in line:
                match = re.search(r'charset (?P<name>[a-zA-Z0-9]+) = (?P<regions>[0-9 -]+);', line)
                regions = []
                for region in match['regions'].split():
                    start, stop = region.split('-')
                    regions.append((int(start)-1, int(stop)))
                transform, start0 = {}, 0
                for start, stop in regions:
                    transform[(start, stop)] = (start0, stop - start + start0)
                    start0 += stop - start
                partition = partitions[partition_id]
                partition.update({'regions': regions, 'transform': transform})
                partition_id += 1

    # Load rate categories
    # In IQ-TREE, only the shape parameter is fit and the rate parameter beta is set to alpha so the mean of gamma distribution is 1
    # The calculations here directly correspond to equation 10 in Yang. J Mol Evol (1994) 39:306-314.
    # Note the equation has a small typo where the difference in gamma function evaluations should be divided by the probability
    # of that category since technically it is the rate given that category
    for partition in partitions.values():
        pinv, alpha, num_categories = partition['pinv'], partition['alpha'], partition['num_categories']
        igfs = []  # Incomplete gamma function evaluations
        for i in range(num_categories+1):
            x = gamma.ppf(i/num_categories, a=alpha, scale=1/alpha)
            igfs.append(gammainc(alpha+1, alpha*x))
        rates = [(0, pinv)]
        for i in range(num_categories):
            rate = num_categories/(1-pinv) * (igfs[i+1] - igfs[i])
            rates.append((partition['speed'] * rate, (1-pinv)/num_categories))
        partition['rates'] = rates

    # Calculate likelihoods
    msa = read_fasta(f'../asr_aa/out/{OGid}.mfa')
    for partition in partitions.values():
        # Unpack partition parameters and partition MSA
        matrix, freqs = models[partition['model']]
        rates = partition['rates']
        partition_msa = []
        for header, seq in msa:
            partition_seq = ''.join([seq[start:stop] for start, stop in partition['regions']])
            partition_msa.append((header, partition_seq))

        # Convert to vectors at tips of tree
        tips = {tip.name: tip for tip in tree.tips()}
        for header, seq in partition_msa:
            tip = tips[header[1:5]]
            conditional = np.zeros((len(alphabet), len(seq)))
            for j, sym in enumerate(seq):
                if sym in sym2idx:
                    i = sym2idx[sym]
                    conditional[i, j] = 1
                else:  # Use uniform distribution for ambiguous symbols
                    conditional[:, j] = 1 / len(alphabet)
            tip.conditional = conditional

        # Calculate likelihoods
        likelihoods = []

        # Get likelihood for invariant category
        # (Background probability for symbol if invariant; otherwise 0)
        likelihood = np.zeros((len(alphabet), len(partition_msa[0][1])))
        for j in range(len(partition_msa[0][1])):
            is_invariant = True
            sym0 = msa[0][1][j]
            if sym0 not in sym2idx:  # Gaps or ambiguous characters are not invariant
                is_invariant = False
            else:
                for i in range(1, len(partition_msa)):
                    sym = partition_msa[i][1][j]
                    if sym != sym0:
                        is_invariant = False
                        break
            if is_invariant:
                idx = sym2idx[sym0]
                likelihood[idx, j] = freqs[idx]
        likelihoods.append(likelihood * rates[0][1])  # Multiply by prior for category

        # Get likelihoods for rate categories
        for rate, prior in rates:
            s, conditional = get_conditional(tree, rate * matrix)
            l = np.expand_dims(freqs, -1) * conditional
            likelihoods.append(np.exp(s) * l * prior)

        likelihoods = np.stack(likelihoods)
        partition['likelihoods'] = likelihoods / likelihoods.sum(axis=(0, 1))

    # Concatenate partition likelihoods
    regions = []
    for partition_id, partition in partitions.items():
        regions.extend([(partition_id, start, stop) for start, stop in partition['regions']])
    regions = sorted(regions, key=lambda x: x[1])

    concatenate = []
    for partition_id, start0, stop0 in regions:
        start, stop = partitions[partition_id]['transform'][(start0, stop0)]
        likelihoods = partitions[partition_id]['likelihoods']
        concatenate.append(likelihoods[:, :, start:stop])
    concatenate = np.concatenate(concatenate, axis=2)

    np.save(f'out/{OGid}_aa.npy', concatenate)

"""
NOTES
The rates for each submodel are calculated manually because some non-invariant submodels are rounded to 0 in IQ-TREE's
output. This results in submodels with zero probability, which introduces problems when normalizing. I felt it was
important to preserve the original model structure when calculating the ASRs, so I decided against merging these
submodels with the invariant submodel.

DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../asr_aa/asr_aa.py
    ../asr_aa/out/*.iqtree
    ../asr_aa/out/*.mfa
    ../asr_aa/out/*.nex
../config/50red_D.txt
../config/WAG.txt
"""