#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
import pandas as pd
import scipy.stats as stats
from scipy.stats import hypergeom
import numpy as np
from statsmodels.stats.multitest import multipletests

all_bed = sys.argv[1]
clade_file = sys.argv[2]
pfam_anno = sys.argv[3]
output_file = sys.argv[4]

def get_according_intersection(left, right, col, out, how): # Vlookup function
    data1 = pd.read_csv(left, sep="\t", low_memory=False)
    data2 = pd.read_csv(right, sep = "\t", low_memory=False)
    data2[col] = data2[col].astype(str)
    data1[col] = data1[col].astype(str)
    df = pd.merge(data1, data2, on=col, how=how)
    df.to_csv(out, index=False, header = True, sep = "\t")

content = []
with open(pfam_anno, 'r') as file:
    for line in file:
        if line.startswith("#"):
            pass
        else:
            c = line.strip().split()
            content.append(c)

df = pd.DataFrame(content) # to form the dataframe
# print(df.head(10))
# print(f"Number of columns: {len(df.columns)}")

# for GN annotation
df_sorted = df.sort_values(by=df.columns[6], ascending=True)  # sort by E-value
df_deduped = df_sorted.drop_duplicates(subset=df.columns[3], keep='first') # uniq the gene_name
columns_to_export = df_deduped.iloc[:, [0, 1, 3, 6, 7, 8]]


columns_to_export.columns = ['Target_name', 'accession', 'geneID', 'E-value', 'score', 'bias']
allbed_column_names = ['Chr', 'Start', 'End', 'Strand', 'geneID']
allbed = pd.read_csv(all_bed, sep='\t',low_memory=False, names=allbed_column_names)
gene_anno = pd.merge(columns_to_export, allbed, on="geneID", how='inner')
columns_to_extract = ['Target_name', 'accession', 'geneID', 'E-value', 'score', 'bias']
gene_anno_to_export = gene_anno[columns_to_extract]

gene_anno_to_export.to_csv("hmmout_anno.txt", sep='\t', header=True, index=False)

# delete the annotation of isoforms
value_counts = columns_to_export.iloc[:, 0].value_counts() # all pfam list
value_counts_dict = value_counts.to_dict() # equal to file_n
N = sum(value_counts_dict.values()) # sum



with open(clade_file, 'r') as r, open("Clade_tmp", 'w') as w:
    w.write("Id\tChr\tStart\tEnd\tStrand\tgeneID\n")
    for line in r:
        w.write(line)

get_according_intersection("Clade_tmp", "hmmout_anno.txt", "geneID", "all.bed.anno.txt", "left")

clade_df = pd.read_csv("all.bed.anno.txt", sep='\t', header=None, skiprows=1) # 
clade_value_counts = clade_df.iloc[:, 6].value_counts() # 第七列
clade_value_counts_dict = clade_value_counts.to_dict() # k字典
K = sum(clade_value_counts_dict.values()) # 被选择K


p_list = []
p_dic = {}
p_dic["Pfam"] = []
p_dic["Odds_ratio"] = []
p_dic["Pvalue"] = []
for key in clade_value_counts_dict:
    k = int(clade_value_counts_dict[key])
    Seleted_list = [k, K-k]
    M = int(value_counts_dict[key])
    Unselected_list = [M, N-M]
    odds_ratio, p_value = stats.fisher_exact([Seleted_list, Unselected_list], alternative='greater')
    p_list.append(p_value)
    p_dic["Pfam"].append(key)
    p_dic["Odds_ratio"].append(odds_ratio)
    p_dic["Pvalue"].append(p_value)
df = pd.DataFrame(p_dic)
alpha = 0.05
corrected_p_values = multipletests(p_list, alpha=alpha, method='bonferroni')[1]
df['adjusted_p_values'] = corrected_p_values
df.to_csv(output_file, index=False, sep='\t')
