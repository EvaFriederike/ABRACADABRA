#!/usr/bin/env python3

import sys
import pandas as pd

input = sys.argv[1]
usher = pd.read_csv(input,sep=',', index_col=0)
usher["FP"] = [0]*len(usher)
usher.to_csv("mod_"+input, sep=',')