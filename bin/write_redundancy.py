#!/usr/bin/env python3

import sys
import os
import utils
import pandas
import json
from collections import defaultdict

fasta = sys.argv[1]
dupl = sys.argv[2]
out = sys.argv[3]
o = open(out, 'w')

D = {}
with open(dupl,'r') as f:
    lines = f.readlines()
    for line in lines:
        ids = [elem.strip(' \n') for elem in line.split('\t')[1].split(',')]
        R = ids[0]
        for header in ids[1:]:
            D[header] = R

for header, seq in utils.parse_fasta(fasta):
    if header in D.keys():
        o.write(f">{D[header]}\n{seq}\n")
o.close()