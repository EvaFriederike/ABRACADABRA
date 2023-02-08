#!/usr/bin/env python3


import sys
import pandas as pd

df_vcf = pd.read_csv(sys.argv[1], sep='\t')

df_vcf.loc[(21563<=df_vcf.POS) & (df_vcf.POS<=25384)] # only consider mutations of the spike region 
df_vcf = df_vcf.loc[(df_vcf["REF"].str.len()==1)&(df_vcf["ALT"].str.len()==1)] # ignore indels

df_vcf.to_csv(sys.argv[1], sep='\t', index=False)