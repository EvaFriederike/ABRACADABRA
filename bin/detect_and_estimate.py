#!/usr/bin/env python3

import sys
import pandas as pd


# read input
vcf = sys.argv[1]
usher = sys.argv[2]
unknown = sys.argv[3]
output = sys.argv[4]

# prepare vcf df structure 
df_vcf = pd.read_csv(vcf, sep='\t')

# if variant_caller == 'medaka':
#     df_vcf.columns = df_vcf.columns.str.upper()
#     df_vcf["POS"] = df_vcf["POS"]+1 # medaka tools vcf2tsv converts 1-based variant positions from vcf file to 0-based coordiantes in tsv file
#     tmp = df_vcf['DPS'].str.split(',')
#     dp = []
#     for r in tmp:
#         dp.append(max(r, key=lambda x: float(x)))
#     df_vcf['DP'] = dp
# if variant_caller == 'ivar':
    #df_vcf.rename(columns={'ALT_DP': 'DP'}, inplace=True)

# df_vcf.loc[(21563<=df_vcf.POS) & (df_vcf.POS<=25384)] # only consider mutations of the spike region 
# df_vcf = df_vcf.loc[(df_vcf["REF"].str.len()==1)&(df_vcf["ALT"].str.len()==1)] # ignore indels
df_usher = pd.read_csv(usher, sep=',', index_col=0)
df_usher.index.name = 'lineage'
df_vcf["var_key"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]
df_vcf['DP'] = df_vcf['ALT_DP'].astype(float)
cluster_muts = df_vcf[["var_key", 'DP']]
cluster_muts = cluster_muts.drop_duplicates(subset='var_key')
cluster_muts.set_index('var_key', inplace=True)
cluster_muts = cluster_muts.transpose()

# prepare usher df structure: consider only mutations that were called 
df_usher_sub = df_usher[["FP"]+[col for col in cluster_muts.columns if col in df_usher.columns]]


## remove variants that exclusively were called for "false positive" lineages
# fp_df = df_usher_sub.loc[df_usher_sub["FP"]==1]
# fp_colsums = fp_df.sum()
# fp_vars = fp_colsums.loc[fp_colsums==1].index.values
# cluster_muts = cluster_muts[[col for col in cluster_muts.columns if col not in fp_vars]]
# print(f"Removing {len(fp_vars)} false positive variants")
# print(f"Remaining variants after removing false positives: {len(cluster_muts.columns)}")
# print(cluster_muts)
# if len(cluster_muts.columns) == 0:
#     sys.exit()
#df_usher_sub = df_usher[["FP"]+[col for col in cluster_muts.columns if col in df_usher.columns]]


# remove false positive lineages
if df_usher_sub['FP'].sum() > 0:
    print(f"Removing false positive lineages")
df_usher_sub = df_usher_sub.loc[df_usher_sub["FP"]==0]
#print(f"Remaining reference lineages after removing potential false positives: {len(df_usher_sub)}")
# add a row for "unknown" lineage
df_usher_sub = df_usher_sub.drop("FP", axis=1)
if unknown == 'true':
    unknown_row = pd.DataFrame(data=[[0]*df_usher_sub.shape[1]], index=["unknown"], columns=df_usher_sub.columns)
    df_usher_sub = pd.concat([df_usher_sub, unknown_row])
    N = cluster_muts.sum(axis=1).item()
else:
    N = df_usher_sub.shape[1]

# for every called mutation, insert read count into corresponding usher columns and normalize by the number of lineages sharing the mutation
# if mutation not in usher table, allocate to "unknown" lineage
print("Compute abundance for each lineage for each mutation by normalizing for collective mutational depth and number of lineages sharing a mutation")
for mut in cluster_muts.columns:
    DP = cluster_muts[mut]
    if isinstance(DP, pd.DataFrame):
        DP = DP.iloc[0,0]
    else:
        DP = DP.item()
    if mut in df_usher_sub.columns:
        colsum = df_usher_sub[mut].sum().item()     
        if colsum != 0:
            print(f"> known mut {mut} from lineages:  {df_usher_sub.loc[df_usher_sub[mut]!=0].index.values}")
            df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP/float(colsum))
        elif unknown == 'true':
            print(f"> unknown mutation {mut}: {DP}") 
            col = [0]*(len(df_usher_sub)-1)+[1]
            tmp = pd.DataFrame(data=[[c] for c in col], columns=[mut])
            tmp.index = df_usher_sub.index
            df_usher_sub = pd.concat([df_usher_sub, tmp], axis=1)
            df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)
    elif unknown == 'true':
        print(f"> unknown mutation {mut}: {DP}")
        col = [0]*(len(df_usher_sub)-1)+[1]
        tmp = pd.DataFrame(data=[[c] for c in col], columns=[mut])
        tmp.index = df_usher_sub.index
        df_usher_sub = pd.concat([df_usher_sub, tmp], axis=1)
        df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)

# normalize matrix by the count of all mutations ~ relative abundance per mutation per lineage
df_usher_sub = df_usher_sub/float(N)

# sum relative abundances of lineage-mutations to final lineage abundance estimates
print("Summing mutational abundances per lineage.")
df_usher_sub['estimate'] = df_usher_sub.sum(axis=1)
print(df_usher_sub)
scaled_cluster_abundances = df_usher_sub.loc[df_usher_sub['estimate']!=0]
#scaled_cluster_abundances = scaled_cluster_abundances*float(factor)
print(scaled_cluster_abundances)
scaled_cluster_abundances.to_csv(output, sep='\t')

