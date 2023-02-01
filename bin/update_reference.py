#!/usr/bin/env python3

import sys
import os
import pandas as pd

df = sys.argv[1]
txt = sys.argv[2]
out = sys.argv[3]

ref = pd.read_csv(df, sep=',', index_col=0)
with open(txt,"r") as f:
    fps = [x.strip() for x in f.readlines()]
mod_ref = ref.loc[[line for line in ref.index if line not in fps],:]
print(f"Removed {len(fps)} from {ref.shape[0]} lineages.")
mod_ref.to_csv(out, sep=',')