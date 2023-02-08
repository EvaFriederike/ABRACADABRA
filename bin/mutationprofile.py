#!/usr/bin/env python3

import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plot_df = pd.DataFrame()
path = sys.argv[1]
out = sys.argv[2]

init = True
for tsv in os.listdir(path):
    if tsv.endswith('.tsv'):
        df_vcf = pd.read_csv(f"{path}/{tsv}", sep='\t')
        cluster = tsv.split('.')[0]
        tmp = pd.DataFrame(columns=['Mut',cluster])
        tmp[cluster] = df_vcf['ALT_FREQ']
        tmp["Mut"] = df_vcf["REF"]+df_vcf["POS"].astype(str)+df_vcf["ALT"]    
        if init:
            plot_df = tmp
            init = False
        else:
            plot_df = plot_df.merge(tmp, how='outer')

plot_df.fillna(0.0, inplace=True)
plot_df.set_index('Mut', inplace=True)
plot_df = plot_df.transpose()
plot_df.sort_index(axis=1, inplace=True)

fig = plt.subplots(figsize=(12,10))
sns.heatmap(plot_df, cmap='magma_r')
plt.savefig(f"{out}_mutation_profile.png", format='png')