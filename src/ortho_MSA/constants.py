"""Default constants for trimming functions."""

constants = {
             # CONSERVED REGIONS PARAMETERS
             'CON_FRAC': 0,  # Maximum gap fraction in conserved columns
             'CON_CLOSE': 3,  # Size of closing element
             'CON_MINLEN': 1,  # Minimum number of columns in conserved regions

             # CONSERVED REGIONS TRIMMING PARAMETERS
             'CON_RATE': 1,  # Decay rate of trim signal
             'CON_MINSIG': 1,  # Minimum signal to trim residues (effectively a minimum number of residues)

             # GAP REGIONS PARAMETERS
             'GAP_NUM': 1,  # Maximum non-gap number in gap columns
             'LOCAL_SIGMA': 1,  # Filter size for calculation of local signal
             'NONLOCAL_RATE': 1,  # Decay rate of nonlocal signal

             # GAP REGIONS TRIMMING PARAMETERS
             'GAP_RATE': 1,  # Decay rate of trim signal
             'GAP_MINSIG': 1,  # Minimum signal to trim residues (effectively a minimum number of residues)

             # LOGISTIC CLASSIFIER PARAMETERS
             'W0': 0,  # Intercept
             'W1': 1,  # Length
             'W2': -1,  # Support
             'W3': -1,  # Local bias
             'W4': 1,  # Nonlocal bias
             'THRESHOLD': 0.5,

             }
