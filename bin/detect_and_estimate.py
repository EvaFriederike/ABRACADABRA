#!/usr/bin/env python3

import sys
import pandas as pd


# read input
vcf = sys.argv[1]
usher = sys.argv[2]
#factor = sys.argv[3]
output = sys.argv[3]

df_vcf = pd.read_csv(vcf, sep='\t')
df_vcf["pos"] = df_vcf["pos"]+1 # medaka tools vcf2tsv converts 1-based variant positions from vcf file to 0-based coordiantes in tsv file
df_vcf.loc[(21563<=df_vcf.pos) & (df_vcf.pos<=25384)] # only consider mutations of the spike region 
df_usher = pd.read_csv(usher, sep=',', index_col=0)
df_usher.index.name = 'lineage'

# prepare vcf df structure 
df_vcf["var_key"] = df_vcf["ref"]+df_vcf["pos"].astype(str)+df_vcf["alt"]
tmp = df_vcf['DPS'].str.split(',')
dp = []
for r in tmp:
    dp.append(max(r, key=lambda x: float(x)))
df_vcf['DP'] = dp
df_vcf['DP'] = df_vcf['DP'].astype(float)
cluster_muts = df_vcf[["var_key", 'DP']]
cluster_muts.set_index('var_key', inplace=True)
cluster_muts = cluster_muts.transpose()

# prepare usher df structure: consider only mutations that were called + add a row for "unknown" lineage
df_usher_sub = df_usher[[col for col in cluster_muts.columns if col in df_usher.columns]]
unknown_row = pd.DataFrame(data=[[0]*df_usher_sub.shape[1]], index=["unknown"], columns=df_usher_sub.columns)
df_usher_sub = pd.concat([df_usher_sub, unknown_row])

# for every called mutation, insert read count into corresponding usher columns and normalize by the number of lineages sharing the mutation
# if mutation not in usher table, allocate to "unknown" lineage
print("Compute abundance for each lineage for each mutation by normalizing for collective mutational depth and number of lineages sharing a mutation")
N = cluster_muts.sum(axis=1).item()
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
            #print(df_usher_sub.loc[df_usher_sub[mut]!=0])
        else:
            print(f"> unknown mutation {mut}: {DP}")
            df_usher_sub[mut] = [0]*(len(df_usher_sub)-1)+[1]
            df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)
            print(df_usher_sub[mut])
    else:
        print(f"> unknown mutation {mut}: {DP}")
        df_usher_sub[mut] = [0]*(len(df_usher_sub)-1)+[1]
        df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)
        print(df_usher_sub[mut])
# normalize matrix by the count of all mutations ~ relative abundance per mutation per lineage
df_usher_sub = df_usher_sub/float(N)

# sum relative abundances of lineage-mutations to final lineage abundance estimates
print("Summing mutational abundances per lineage.")
df_usher_sub['estimate'] = df_usher_sub.sum(axis=1)
scaled_cluster_abundances = df_usher_sub.loc[df_usher_sub['estimate']!=0]
#scaled_cluster_abundances = scaled_cluster_abundances*float(factor)
scaled_cluster_abundances.to_csv(output, sep='\t')

