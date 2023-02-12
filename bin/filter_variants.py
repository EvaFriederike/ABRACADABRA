#!/usr/bin/env python3


import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import pysam
from Bio import SeqIO



def ecdf(values):
    x, counts = np.unique(values, return_counts=True)
    cumsum = np.cumsum(counts)
    return x,cumsum/cumsum[-1]


def filter_variants(known_df, unknown_df):
    """
    Filter unknown variants and keep only those where the AF is larger or equal to
    the AF of the known distributio while the ECDF value if smaller or equal to the corresponding
    ECDF values of known variants
    """
    known_af, known_prob = ecdf(known_df["ALT_FREQ"])
    unknown_af, unknown_prob = ecdf(unknown_df["ALT_FREQ"])
    known_data = list(zip(known_af, known_prob))
    unknown_data = list(zip(unknown_af, unknown_prob))
    af_cutoff = known_data[::-1][0][0]
    for tup in known_data[::-1]:
        ovlp = [x[0] for x in unknown_data if x[0]>=tup[0] and x[1]<=tup[1]]
        if len(ovlp) != 0:
            af_cutoff = min(ovlp)
        else:
            break
 
    unknown_df = unknown_df.loc[unknown_df["ALT_FREQ"]>=af_cutoff]

    return pd.concat([known_df,unknown_df])



vcf = sys.argv[1]
usher = sys.argv[2]
unknown = sys.argv[3]
bam = sys.argv[4]
ref = sys.argv[5]
df_vcf = pd.read_csv(vcf, sep='\t')
df_usher = pd.read_csv(usher, sep=',', index_col=0)

# prepare vcf df structure: filter for passing variants and SNPs
df_vcf = df_vcf.loc[df_vcf['PASS']==True]
df_vcf = df_vcf.loc[(df_vcf["REF"].str.len()==1)&(df_vcf["ALT"].str.len()==1)] # ignore indels
df_vcf["var_key"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]
df_vcf = df_vcf.drop_duplicates(subset='var_key')

# prepare usher df structure: consider only mutations that were called 
df_usher_sub = df_usher[[col for col in df_vcf["var_key"].values if col in df_usher.columns]]
df_vcf["Group"] = df_vcf['var_key'].isin(df_usher_sub.columns)
df_vcf["Group"].replace(True, "known", inplace=True)
df_vcf["Group"].replace(False, "unknown", inplace=True)
# Filter unknown variants for potential false positives
# if unknown mode is on 
if unknown == 'true':
    # if there are unknown mutations
    if df_usher_sub.shape[1]!= 0 and len(df_vcf) > df_usher_sub.shape[1]:  
        # if there are more known than unknown mutations:
        #if len(df_vcf[df_vcf["Group"]=="known"]) > len(df_vcf[df_vcf["Group"]=="unknown"]):
        plt.figure(figsize=(10,8))
        sns.histplot(x='ALT_FREQ', data=df_vcf, hue='Group', bins=len(df_vcf), stat='density', element='step', fill=False, cumulative=True, common_norm=False)
        plt.savefig("AF_histogram.png")
        true_vars = filter_variants(df_vcf[df_vcf["Group"]=="known"], df_vcf[df_vcf["Group"]=="unknown"])
        before = len(df_vcf)
        after = len(true_vars)
        print(f"Filtered {before-after} variants as potential false positive calls")

        # get variants where reads that contain them need to be filtered out of the bamfile
        # remove reads only if they carry only an obsolete var 
        if before-after != 0:
            false_vars = df_vcf[~df_vcf["var_key"].isin(true_vars["var_key"])]
            keep = []
            toss = []
            genome = SeqIO.to_dict(SeqIO.parse(open(ref,"r"),"fasta"))
            genome = str(genome['NC_045512.2'].seq)
            bamfile = pysam.AlignmentFile(bam, "rb")  
            for pileupcolumn in bamfile.pileup(reference='NC_045512.2'):
                refpos = pileupcolumn.reference_pos
                #false_vars.POS.values:
                if refpos+1 in df_vcf.POS.values: 
                    #print(f"Checking for reads containing false positive filtered variants at reference position: {refpos}={genome[refpos]}")
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            qpos = pileupread.query_position
                            var = genome[refpos]+str(refpos)+pileupread.alignment.query_sequence[qpos]
                            if var in false_vars["var_key"].values:
                                toss.append(pileupread.alignment.qname)
                            elif var in true_vars["var_key"].values:
                                keep.append(pileupread.alignment.qname)
           
            read_names = [r for r in toss if r not in keep]
            with open("filtered-reads.txt","w") as f:
                for rname in read_names:
                    f.write(rname+'\n')
            df_vcf = true_vars
else:
    df_vcf = df_vcf[df_vcf["Group"]=="known"]

for c in ['var_key','Group']:        
    if c in df_vcf.columns:
        df_vcf = df_vcf.drop(c, axis=1)
    
df_vcf.to_csv("filtered_"+vcf, sep='\t', index=False)