#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import regex as re
import utils
import json
from collections import defaultdict
import pandas as pd
import numpy as np

# read input data
sequenceFile = sys.argv[1]
primerFile = sys.argv[2]
readcount = sys.argv[3]
mismatch_cutoff = sys.argv[4]
df = pd.read_csv(primerFile, sep='\t')


tmp = defaultdict(list)
none = {}
# compare every read against every primer seq allowing 1-2 mismatches
QC= 0
sorted_reads = 0
split_reads = 0
subreads = 0
for header, sequence in utils.parse_fasta(sequenceFile):
    QC += 1
    match_count = 0
    matching_primer = []
    for i,r in df.iterrows():
        m = re.findall('('+r.seq+'){s<='+mismatch_cutoff+'}', sequence)
        if len(m) != 0:
            match_count += 1
            matching_primer.append(r['name'])
    # reads matching exactly one primer are assigned to the corresponding primer set
    if match_count == 1:
        sorted_reads += 1
        tmp[matching_primer[0]].append([header, sequence])
    # reads matching no primer are stored
    elif match_count == 0:
        none[header] = sequence
    # reads matching >= 2 primers are split into sub-reads based on the starting position of the contained reads 
    elif match_count >=2:
        split_reads += 1
        start_pos = {}
        c = 0
        # get starting positions of every matched primer in the read
        for primer in matching_primer:
            primer_seq = df.loc[df['name']==primer].seq.item()
            s = re.search('('+primer_seq+'){s<='+mismatch_cutoff+'}', sequence).start()
            start_pos[primer] = s
        # catch first subread (this includes matching primer AND some remaining adapter seq)
        start_pos_sorted = sorted(start_pos.items(), key= lambda x: x[1])
        init = start_pos_sorted[1]
        subsequence = sequence[:init[1]]
        subreads += 1
        tmp[init[0]].append([f"{header}_sub{c}", subsequence])
        # print(f"- start {start_pos_sorted[0]} to {start_pos_sorted[1]}: len={len(subsequence)}")
        start_pos_sorted = start_pos_sorted[1:]
        c += 1
        # split remaining subreads
        for k,v in start_pos_sorted[::-1]:
            subreads += 1
            tmp[k].append([f"{header}_sub{c}", sequence[v:]])
            # print(f"- start ({k},{v}): len={len(sequence[v:])}")
            subsequence = subsequence[:s]
            c += 1


log = open('primer_sort.log', 'w')
log.write(f"{QC}/{readcount} ({round(QC/float(readcount)*100,2)}%) sample reads passed N-based QC.\n")
log.write(f"{sorted_reads+split_reads} ({round((sorted_reads+split_reads)/float(readcount)*100,2)}%) of all input reads could be sorted into an amplicon bin. {len(none)} reads could not be assigned.\n")
log.write(f"{split_reads} reads matched more than one amplicon and had to be split, resulting in a total of {sorted_reads+subreads} amplicon reads.\n")
log.write(f"--> {sorted_reads+subreads} reads were sorted into {len(set([ k.split('_')[1] for k in tmp.keys()]))} of {int(len(df)/2)} amplicons. <--\n")
log.write(f"On average, an amplicon contains {int(np.mean(np.array([len(v) for v in tmp.values()])))} reads:\n")
sorted_tmp = sorted(tmp.items(), key= lambda x: len(x[1]))
log.write(f"Amplicon {sorted_tmp[0][0].split('_')[1]} has the smallest number of reads with {len(sorted_tmp[0][1])}.\nAmplicon {sorted_tmp[-1][0].split('_')[1]} has the largest number of reads with {len(sorted_tmp[-1][1])}.\n")
log.close()

# write output fasta files   
for k,v in tmp.items():
    f = open(f"{k}.fasta", 'w')
    for r in v:
        f.write(f">{r[0]}\n{r[1]}\n")
    f.close()

f = open(f"no_primer_match.fasta", 'w')
for k,v in none.items():  
    f.write(f">{k}\n{v}\n")
f.close()
