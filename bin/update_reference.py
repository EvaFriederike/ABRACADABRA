#!/usr/bin/env python3

import sys
import os
import pandas as pd

df = sys.argv[1]
out = sys.argv[2]

ref = pd.read_csv(df, sep=',', index_col=0)

ref["FP"] = [0]*len(ref)

ref.to_csv(out, sep=',')