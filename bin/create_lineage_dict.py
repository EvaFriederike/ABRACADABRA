#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import utils
# import json
import pandas as pd

sequenceFile = sys.argv[1]
lineage = sys.argv[2]

D = {}

for header, sequence in utils.parse_fasta(sequenceFile):
   D[header] = lineage


lineagedict = pd.DataFrame(data=D.items(), columns=["header","lineage"])
lineagedict.to_csv("lineageDict.csv", sep=',', index=False)