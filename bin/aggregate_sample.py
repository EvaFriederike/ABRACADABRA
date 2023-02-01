#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np


# read input
tsv = sys.argv[1]
T = float(sys.argv[2])
output = sys.argv[3]

df = pd.read_csv(tsv, sep='\t', index_col=0)

relabs = []
lines = []
for lineage in df.index.unique():
    relab = df.loc[lineage].sum().item()
    lines.append(lineage)
    relabs.append(relab)

final_abundances = pd.DataFrame(data=[relabs], columns=lines)
final_abundances = final_abundances.transpose()
final_abundances.columns = ['relative abundance']
fps = final_abundances.loc[final_abundances["relative abundance"]<T].index.values
if len(fps) > 0:
    print(f"{len(fps)} lineages were estimated with an abundance below the false positive threshold!\nRemove probable false positives from reference and recalculate:")
    with open("FP_lineages.txt","w") as f:
        for l in fps:
            print(f"{l}\n")
            f.write(l+'\n')
else:
    print(f"The final lineage abundances sum to {round(final_abundances.sum().item()*100,2)}% ")
final_abundances.to_csv(output, sep='\t')