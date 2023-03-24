#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Source: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
# "Clustering of viral genomes based on different algorithms and metrices"
# (Author: kevin.lamkiewicz@uni-jena.de)

import sys

import utils

sequenceFile = sys.argv[1]
outputFile = sequenceFile.replace("positive","negative")
print(sequenceFile)
with open(outputFile, 'w') as outputStream:
  for header, sequence in utils.parse_fasta(sequenceFile):
    print(header,sequence)
    outputStream.write(f">{header}\n{utils.reverseComplement(sequence)}\n")
