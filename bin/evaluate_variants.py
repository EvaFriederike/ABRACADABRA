#!/usr/bin/env pyhton3

import sys
import pandas as pd

var_df = pd.read_csv(sys.argv[1], sep='\t')
usher = pd.read_csv(sys.argv[2], index_col=0)
if sys.argv[3] != "dummy":
    sample = pd.read_csv(sys.argv[3])
    lineages = sample.lineage.unique()
    usher = usher.loc[lineages,:]

var_df['VAR'] = var_df['REF']+var_df['POS'].astype(str)+var_df['ALT']
var_df.drop_duplicates('VAR', inplace=True)
usher_sub = usher[usher["FP"]==0]

TP = [col for col in var_df.VAR if col in usher_sub.columns.values]
known_vars = len(TP)/float(len(var_df))
unknown_vars = (len(var_df)-len(TP))/float(len(var_df))

with open('variant_evaluation.log','w') as f:
    f.write(f"Proportion of known called variants: {known_vars}\n")
    f.write(f"Proportion of unknown called variants: {unknown_vars}\n")
    if unknown_vars != 0:
        f.write(f"unknown vars:\n")
        for v in [x for x in var_df["VAR"] if x not in TP]:
            f.write(f"\t{v}\n")