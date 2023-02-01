#!/usr/bin/env python3

import sys
import pandas as pd


# read input
tsv = sys.argv[1]
factor = float(sys.argv[2])
output = sys.argv[3]

df = pd.read_csv(tsv, sep='\t', index_col=0)
scaled_df = df*float(factor)
scaled_df.to_csv(output, sep='\t')

