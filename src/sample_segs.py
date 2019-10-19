"""Create samples sequences at different cutoffs."""

import pandas as pd
from sys import argv

# Input variables
path = argv[1]
type_name = argv[2]  # Name of column denoting segment type

# Constants
n = 8500

segs = pd.read_csv(path, sep='\t', keep_default_na=False)
segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
for i in [2 ** x for x in range(6)]:
    # Read data and split segments
    segs_lt = segs[segs['seq'].map(lambda x: len(x) >= i)]
    T_segs = segs_lt[segs_lt[type_name]]
    F_segs = segs_lt[~segs_lt[type_name]]

    # Sample segments
    T_sample = T_segs.sample(n).sort_values('seg_id')
    F_sample = F_segs.sample(n).sort_values('seg_id')

    # Merge samples and save
    sample = pd.concat([T_sample, F_sample])
    sample.to_csv(f'segments_{i}.tsv', sep='\t', index=False)
