#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import utils
import regex


sequenceFile = sys.argv[1]
lprimer = sys.argv[2]
rprimer = sys.argv[3]
m = sys.argv[4]
f = open('split_reads.fasta','w')

lregex = regex.compile(r"("+lprimer+"){s<="+m+"}")
rregex = r"("+rprimer+"){s<="+m+"}"

for header, sequence in utils.parse_fasta(sequenceFile):
    # loop over finditer matches for lprimer
    # collect leftmost span index and if presentn_substitutions
    # sort matches by left most span and secondary ba n_substitutions descending
    # pick first item span as cut position
    lmatch_pos = []
    for match in regex.finditer(lregex, sequence):
        lmatch_pos.append(match.span()[0])
    lmatch_pos.sort()
    rmatch_pos = []
    for match in regex.finditer(rregex, sequence):
        rmatch_pos.append(match.span()[0])
    rmatch_pos.sort()

    if len(lmatch_pos) != 0 and len(rmatch_pos) != 0:
        l_start = lmatch_pos[0]
        r_end = rmatch_pos[-1]
        read = sequence[l_start:r_end+1]
    elif len(lmatch_pos) == 0:
        if len(rmatch_pos) == 0:
            #print(f"Weird.....read {header} doesn't match any of the amplicon primers it matches before using seqkit grep?")
            read = sequence
        else:
            r_end = rmatch_pos[-1]
            read = sequence[:r_end+1]
    else:
        l_start = lmatch_pos[0]
        read = sequence[l_start:]
    f.write(f">{header}\n{read}\n")
    







