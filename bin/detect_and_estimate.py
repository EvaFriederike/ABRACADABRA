#!/usr/bin/env python3

import sys
import pandas as pd

# read input
vcf = sys.argv[1]
usher = sys.argv[2]
unknown = sys.argv[3]
output = sys.argv[4]

df_usher = pd.read_csv(usher, sep=',', index_col=0)

# prepare vcf df structure 
df_vcf = pd.read_csv(vcf, sep='\t')
df_vcf["var_key"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]
df_vcf['DP'] = df_vcf['ALT_DP'].astype(float)
cluster_muts = df_vcf[["var_key", 'DP']]
cluster_muts.set_index('var_key', inplace=True)
cluster_muts = cluster_muts.transpose()

if df_usher["FP"].sum() != 0:
    print("Removing false positive lineages")
df_usher_sub = df_usher.loc[df_usher["FP"]==0]
df_usher_sub = df_usher[[col for col in cluster_muts.columns if col in df_usher_sub.columns]]

# if unknown mode is on 
if unknown == 'true':
    # add a row for "unknown" lineage
    unknown_row = pd.DataFrame(data=[[0]*df_usher_sub.shape[1]], index=["unknown"], columns=df_usher_sub.columns)
    df_usher_sub = pd.concat([df_usher_sub, unknown_row])
else:
    cluster_muts = cluster_muts[[col for col in cluster_muts.columns if col in df_usher_sub.columns]]

N = cluster_muts.sum(axis=1).item()
print(f"N={N}, n_cluster_muts={cluster_muts.shape[1]}, n_shared_muts={df_usher_sub.shape[1]}")
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
            print(f"> known mutation but no matching lineage in UShER dataframe {mut}: {DP}") 
            col = [0]*(len(df_usher_sub)-1)+[1]
            tmp = pd.DataFrame(data=[[c] for c in col], columns=[mut])
            tmp.index = df_usher_sub.index
            df_usher_sub = pd.concat([df_usher_sub, tmp], axis=1)
            df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)
        else:
            print(f"> known mutation but no matching lineage in UShER dataframe {mut}: {DP}")
            print(f"WARNING: DP is excluded from estimation, but associated reads are still be used for scaling! Final abundances might not sum to 100%. Please check that you remove mutations from UShER dataframe that are not captured by any lineage!")
            N -= DP
    elif unknown == 'true':
        print(f"> unknown mutation {mut}: {DP}")
        col = [0]*(len(df_usher_sub)-1)+[1]
        tmp = pd.DataFrame(data=[[c] for c in col], columns=[mut])
        tmp.index = df_usher_sub.index
        df_usher_sub = pd.concat([df_usher_sub, tmp], axis=1)
        df_usher_sub[mut] = df_usher_sub[mut].where(df_usher_sub[mut]==0, DP)
    else:
        print(f"WARNING: unknown mutation {mut} is still in the mix despite running in known-only mode. DP is used for estimation")
# normalize matrix by the count of all mutations ~ relative abundance per mutation per lineage
df_usher_sub = df_usher_sub/float(N)

# sum relative abundances of lineage-mutations to final lineage abundance estimates
print("Summing mutational abundances per lineage.")
df_usher_sub['estimate'] = df_usher_sub.sum(axis=1)
scaled_cluster_abundances = df_usher_sub.loc[df_usher_sub['estimate']!=0]

scaled_cluster_abundances.to_csv(output, sep='\t')