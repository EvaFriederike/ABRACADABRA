#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import utils
import json

sequenceFile = sys.argv[1]
primerDict = json.load(open(sys.argv[2], 'r'))

L = open("LEFT.fasta", "w")
R = open("RIGHT.fasta", "w")

for header, sequence in utils.parse_fasta(sequenceFile):
    if primerDict[header] == 'LEFT':
        L.write(f"\n>{header}\n{sequence}")
    else:
        R.write(f"\n>{header}\n{sequence}")