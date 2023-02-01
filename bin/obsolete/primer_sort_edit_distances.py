#!/usr/bin/env python3

import sys
import utils
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict


def calc_edit_distance(seq, query):
    """
    """
    ed_list = []
    for i in range(len(seq)-len(query)):
        ed = 0
        c = 0
        for j in range(i, i+len(query)):
            if seq[j] != query[c]:
                ed += 1
            c += 1
        ed_list.append(ed/float(len(query)))
    min_ed = min(ed_list, default='EMPTY')

    return min_ed

# read primer df and fasta with no primer matches
primer = sys.argv[1]
fasta = sys.argv[2]
df = pd.read_csv(primer, sep='\t')

# for s in fasta: 
#   for p in primer:
#       calc min edit distance between p and s
#   take min edit distance between s and p and memorize p
ed_values = []
primer_list = []
for header, sequence in utils.parse_fasta(fasta):
    ed_dict = defaultdict(list)
    for i,r in df.iterrows():
        ed = calc_edit_distance(sequence, r.seq)
        ed_dict[r['name']].append(ed)
    min_ed = min([x for k,v in ed_dict.items() for x in v])
    min_ed_primer = [k for k,v in ed_dict.items() if min_ed in v]
    print(f"ed({header},{min_ed_primer[:2]},...) = {min_ed}")
    
    ed_values.append(min_ed)
    primer_list.extend(min_ed_primer)
    
    

# plot histogram over primers that show min edit distance against a read
# plot histogram for the observed min edit distance values
plt.hist(ed_value)
plt.savefig('min-edit-distance.png')

plt.hist(min_ed_primer)
plt.savefig('primers-with-min-edit-distance.png')