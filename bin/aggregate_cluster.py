#!/usr/bin/env python3

import sys
import os
import pandas as pd


# read input
path = sys.argv[1]
output = sys.argv[2]

files = [file for file in os.listdir(path)]
init_df = pd.read_csv(path+files[0], sep='\t', index_col=0)
init_df = init_df['estimate']

for file in files[1:]:
    df = pd.read_csv(path+file, sep='\t', index_col=0)
    init_df = pd.concat([init_df, df['estimate']])

lines = []
relabs = []
for lineage in init_df.index.unique():
    relab = init_df[lineage].sum().item()
    relabs.append(relab)
    lines.append(lineage)
    print(f"summing {lineage} across clusters: {relab}")

amplicon_abundances = pd.DataFrame(data=[relabs], columns=lines)
amplicon_abundances = amplicon_abundances.transpose()
amplicon_abundances.columns = ['estimate']
amplicon_abundances.to_csv(output, sep='\t')