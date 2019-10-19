"""Create samples sequences at different cutoffs."""

import pandas as pd
from sys import argv

# Input variables
path = argv[1]

# Constants
n = 8500

segs = pd.read_csv(path, sep='\t', keep_default_na=False)
segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
for i in [2 ** x for x in range(6)]:
    segs_lt = segs[segs['seq'].map(lambda x: len(x) >= i)]
    sample = segs_lt.sample(n)
    sample.to_csv(f'segments_{i}.tsv', sep='\t', index=False)
