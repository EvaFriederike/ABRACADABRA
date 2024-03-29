#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np


# read input
tsv = sys.argv[1]
T = float(sys.argv[2])
ref = sys.argv[3]
unknown = sys.argv[4]

usher = pd.read_csv(ref, sep=',', index_col=0)
df = pd.read_csv(tsv, sep='\t', index_col=0)

relabs = []
lines = []
for lineage in df.index.unique():
    relab = df.loc[lineage].sum().item()
    print(f"{lineage}: {relab}")
    lines.append(lineage)
    relabs.append(relab)

final_abundances = pd.DataFrame(data=[relabs], columns=lines)
final_abundances = final_abundances.transpose()
final_abundances.columns = ['relative abundance']
if usher["FP"].sum()==0:
    fps = final_abundances.loc[final_abundances["relative abundance"]<T].index.values
    if 'unknown' in fps:
            fps = fps[fps!='unknown']
    if len(fps) > 0:
        output = "unfiltered_sample_abundances.tsv"
        # update reference for false positives
        usher.loc[fps,['FP']] = 1
        usher.to_csv("barcodes_FP_annotated.csv", sep=',')
        with open('false_positive_detection.log', 'w') as f:
            f.write(f"{len(fps)} lineages were estimated with an abundance below the false positive threshold!\nRemove probable false positives from reference and recalculate:\n")
            for l in fps:
                f.write(f"{l}\n")
    else:
        output = "final_sample_abundances.tsv"
        print(f"The final lineage abundances sum to {round(final_abundances.sum().item()*100,2)}% ")
else:
    output = "final_sample_abundances.tsv"
    print(f"The final lineage abundances sum to {round(final_abundances.sum().item()*100,2)}% ")
    # if unknown == 'true' and final_abundances.loc['unknown','relative abundance'] < T:
    #     final_abundances.drop(index='unknown', inplace=True)

final_abundances.to_csv(output, sep='\t')

