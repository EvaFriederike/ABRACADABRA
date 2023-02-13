#!/usr/bin/env python3


import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import pysam
from Bio import SeqIO

def filter_reads(true_vars, ref, bam):
    """
    Keep only reads that contain true positive variants
    """
    keep = []
    genome = SeqIO.to_dict(SeqIO.parse(open(ref,"r"),"fasta"))
    genome = str(genome['NC_045512.2'].seq)
    bamfile = pysam.AlignmentFile(bam, "rb")  
    for pileupcolumn in bamfile.pileup(reference='NC_045512.2'):
        refpos = pileupcolumn.reference_pos
        if refpos+1 in true_vars.POS.values: 
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    qpos = pileupread.query_position
                    var = genome[refpos]+str(refpos+1)+pileupread.alignment.query_sequence[qpos]
                    if var in true_vars["var_key"].values:
                        keep.append(pileupread.alignment.qname)

    read_names = set(keep)
    with open("filtered-reads.txt","w") as f:
        for rname in read_names:
            f.write(rname+'\n')
            
    return None

############################################### MAIN ##############################################


vcf = sys.argv[1]
usher = sys.argv[2]
unknown = sys.argv[3]
bam = sys.argv[4]
ref = sys.argv[5]
df_vcf = pd.read_csv(vcf, sep='\t')
df_usher = pd.read_csv(usher, sep=',', index_col=0)
df_vcf["var_key"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]
print(f"Before filtering: {len(df_vcf)} variants")

# get all variants from vcf that are associated with false positive lienages in the usher frame
# remove fp lineages and mutations with colusm 0 from usher
# remove selected variants from vcf table
# loop through bam file and identify reads that carry removed mutations
df_usher_sub = df_usher.loc[df_usher["FP"]==0]
# identify mutations that are now without lineage assignment and remove them from the table of called mutations
# since those are now "unknown"
obsolete_vars = [col for col in df_usher_sub.columns if df_usher_sub[col].sum()==0] 
var_pos = df_vcf.POS.values
df_vcf = df_vcf[~df_vcf['var_key'].isin(obsolete_vars)]
print(f"After filtering: {len(df_vcf)} variants")

# df_usher_sub = df_usher_sub[[c for c in df_usher_sub.columns if not c in obsolete_vars]]
# df_usher_out = df_usher_sub.merge(df_usher['FP'], left_index=True, right_index=True)

if len(df_vcf) != 0:
    filter_reads(df_vcf, ref, bam)

df_vcf = df_vcf.drop('var_key', axis=1)
df_vcf.to_csv("filtered_"+vcf,sep='\t')
# df_usher_out.to_csv("mod_"+usher, sep=',')