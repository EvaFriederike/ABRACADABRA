#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Source: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
# "Clustering of viral genomes based on different algorithms and metrices"
# (Author: kevin.lamkiewicz@uni-jena.de)

import sys
import re
from collections import defaultdict

genbankACCRegex = re.compile(r'[^A-Z]*([A-Z]{2}[0-9]{8}|[A-Z]{2}[0-9]{6}|[A-Z][0-9]{5}|NC_[0-9]{6})[^0-9]?')

def reverseComplement(sequence):
  """

  """
  comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
  return(''.join([comp.get(x,'N')  for x in sequence.upper()[::-1]]))


def parse_fasta(filePath):
  """
  """
  with open(filePath, 'r') as inputStream:
    header = ''
    seq = ''

    for line in inputStream:
      if line.startswith(">"):
        if header:
          if seq.count('N') / len(seq) < 0.1:
            yield (header, seq)

        header = line.rstrip("\n ").replace(':','_').replace(' ','_').replace('|','_').lstrip(">").rstrip('_')
        seq = ''
      else:
        seq += line.rstrip("\n").upper().replace('U','T')

    if seq.count('N') / len(seq) < 0.1:
      yield (header, seq)
