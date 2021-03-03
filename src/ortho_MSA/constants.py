"""Default constants for trimming functions."""

from importlib.resources import open_text

matrix = {}
with open_text(__package__, 'BLOSUM62.txt') as file:
    file.readline()  # Skip header
    syms = file.readline().split()
    for i, line in enumerate(file):
        for j, value in enumerate(line.split()[1:]):
            matrix[(syms[i], syms[j])] = int(value)

constants = {
             # CONSERVED REGIONS PARAMETERS
             'CON_FRAC': 0,  # Maximum gap fraction in conserved columns
             'CON_WINDOW': 3,  # Size of closing element
             'CON_MINLEN': 1,  # Minimum number of columns in conserved regions

             # CONSERVED REGIONS TRIMMING PARAMETERS
             'CON_RATE': 1,  # Decay rate of trim signal
             'CON_MINSIG': 1,  # Minimum signal to trim residues (effectively a minimum number of residues)

             # GAP REGIONS PARAMETERS
             'GAP_NUM': 1,  # Maximum non-gap number in gap columns
             'GP_SIGMA': 1,  # Filter size for calculation of gap propensity
             'GD_WINDOW': 1,  # Size of gap diversity window
             'INDEL1_RATE': 1,  # Decay rate of indel signal 1
             'INDEL2_RATE': 1,  # Decay rate of indel signal 2
             'MATRIX': matrix,

             # GAP REGIONS TRIMMING PARAMETERS
             'GAP_RATE': 1,  # Decay rate of trim signal
             'GAP_MINSIG': 1,  # Minimum signal to trim residues (effectively a minimum number of residues)

             # LOGISTIC CLASSIFIER PARAMETERS
             'W0': 0,  # Intercept
             'W1': 1,  # Length
             'W2': -1,  # Support
             'W3': -1,  # Gap propensity
             'W4': 1,  # Gap diversity
             'W5': 1,  # Indel bias 1
             'W6': 1,  # Indel bias 2
             'THRESHOLD': 0.5,
             }
