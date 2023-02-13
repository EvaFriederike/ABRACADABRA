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
    print(f"Final ALT_FREQ cutoff for potential false variants: {af_cutoff}")
    unknown_df = unknown_df.loc[unknown_df["ALT_FREQ"]>=af_cutoff]

    return pd.concat([known_df,unknown_df])


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

# prepare vcf df structure: filter for passing variants and SNPs
before = len(df_vcf)
print(f"Before filtering: {before} variants")
df_vcf = df_vcf.loc[df_vcf['PASS']==True]
df_vcf = df_vcf.loc[(df_vcf["REF"].str.len()==1)&(df_vcf["ALT"].str.len()==1)] # ignore indels
df_vcf["var_key"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]
df_vcf = df_vcf.drop_duplicates(subset='var_key')
after = len(df_vcf)
print(f"Removed {before-after} variants that lied outside the S gene or were indels")

# prepare usher df structure: consider only mutations that were called 
df_usher_sub = df_usher[[col for col in df_vcf["var_key"].values if col in df_usher.columns]]
df_vcf["Group"] = df_vcf['var_key'].isin(df_usher_sub.columns)
df_vcf["Group"].replace(True, "known", inplace=True)
df_vcf["Group"].replace(False, "unknown", inplace=True)

# if unknown mode is on, filter unknown variants for potential false positives
if unknown == 'true':
    # if there are known and unknown mutations
    if df_usher_sub.shape[1]!= 0:
        if len(df_vcf) > df_usher_sub.shape[1]:  
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
            filter_reads(true_vars, ref, bam)
            df_vcf = true_vars
        # if there are no unknown mutations, output df_vcf and select reads carrying known variants
        else:
            df_vcf = df_vcf[df_vcf["Group"]=="known"]
            filter_reads(df_vcf, ref, bam)
    # if there are no known variants, output empty vcf dataframe and no selected reads
    else:
        df_vcf = df_vcf[df_vcf["Group"]=="known"]
    
# in known mode and if there are unknown variants, keep only reads that carry known variants
elif unknown == 'false':
    if df_usher_sub.shape[1]!= 0:
        if len(df_vcf) > df_usher_sub.shape[1]:
            true_vars = df_vcf[df_vcf["Group"]=="known"]
            unknown_vars = df_vcf[df_vcf["Group"]=="unknown"]
            print(f"Removed {len(unknown_vars)} unknown variants in known mode")
            filter_reads(true_vars, ref, bam)
            df_vcf = true_vars
        # if there are no unknown variants,  output df_vcf and select reads
        else:
            df_vcf = df_vcf[df_vcf["Group"]=="known"]
            filter_reads(df_vcf, ref, bam)
    # if there are no known variants, output empty vcf dataframe and select no reads
    else:
        df_vcf = df_vcf[df_vcf["Group"]=="known"]
    


for c in ['var_key','Group']:        
    if c in df_vcf.columns:
        df_vcf = df_vcf.drop(c, axis=1)

print(f"After filtering: {len(df_vcf)} variants")
df_vcf.to_csv("filtered_"+vcf, sep='\t', index=False)