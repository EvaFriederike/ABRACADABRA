#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################
# TODO:
# * reads don't start with primers => barcode tags?
# * reads can contain multiple primers, implement a strategy, e.g. only select a forward
#   read for the amplicon of interest, if it doesn't contain the previous forward primer and 
#   vice versa for the reverse strand reads
# * validate with NanoPlot / filter for length?
######################################

import sys
import re
import utils
import json

sequenceFile = sys.argv[1]
primer_left = sys.argv[2]
primer_right = sys.argv[3]

R = open("right_sorted.fasta", "w")
L = open("left_sorted.fasta", "w")
N = open("not_sorted.fasta", "w")
primerDict = open('primerDict.json', 'w')
tmp = {}

for header, sequence in utils.parse_fasta(sequenceFile):
    if primer_left in sequence:
        if primer_right in sequence:
            M.write(f"\n>{header}\n{sequence}")
        else:
            L.write(f"\n>{header}\n{sequence}")
            tmp[header] = 'LEFT'
    elif primer_right in sequence:
        R.write(f"\n>{header}\n{sequence}")
        tmp[header] = 'RIGHT'
    else:
        N.write(f"\n>{header}\n{sequence}")
    
json.dump(tmp, primerDict)
