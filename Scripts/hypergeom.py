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

def get_according_intersection(left, right, col, out, how): # 根据某一列合并两个文件
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

# 和基因表里面进行取交集后再统计，pfam_anno里有
value_counts = columns_to_export.iloc[:, 0].value_counts() # all pfam list
value_counts_dict = value_counts.to_dict() # equal to file_n
N = sum(value_counts_dict.values()) # 累加



with open(clade_file, 'r') as r, open("Clade_tmp", 'w') as w:
    w.write("Id\tChr\tStart\tEnd\tStrand\tgeneID\n")
    for line in r:
        w.write(line)
# os.system("sed -i '1i\Target_name\taccession\tgeneID\tE-value\tscore\tbias\n' hmmout_anno.txt")
# 单独这个Clade的结果
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

"""
version2
p_list = []
p_dic = {}
for key in clade_value_counts_dict:
    k = int(clade_value_counts_dict[key])
    n = K
    M = int(value_counts_dict[key])
    pvalue = 1-hypergeom.cdf(k, N, M, n)
    p_list.append(pvalue)
    p_dic[key] = pvalue
df = pd.DataFrame(list(p_dic.items()), columns=['Key', 'Value'])
_, adjusted_p_values, _, _ = multipletests(p_list, method='fdr_bh')
df['adjusted_p_values'] = adjusted_p_values
df.to_csv(output_file, index=False, sep='\t')

file_K = sys.argv[1] # 被选择成功
file_n = sys.argv[2] # 总体中的成功
N = int(sys.argv[3]) # 总体大小
output_file = sys.argv[4] # 这里直接生成

version 1
p_list = []
p_dic = {}
with open(file_K, 'r') as r, open(output_file, 'w') as w:
    r.readline()
    num = 0
    for i in r:
        c = i.strip().split("\t")
        pfam_id = c[0]
        k = int(c[1])
        n = int(c[2])
        with open(file_n, 'r') as r2:
            r2.readline()
            for line in r2:
                C = line.strip().split('\t')
                if C[0] == pfam_id:
                    M = int(C[1])
                    break
        # print("1-hypergeom.cdf({}, {}, {}, {})".format(k, N, M, n))
        pvalue = 1-hypergeom.cdf(k, N, M, n)
        p_list.append(pvalue)
        p_dic[pfam_id] = pvalue
    df = pd.DataFrame(list(p_dic.items()), columns=['Key', 'Value'])
    _, adjusted_p_values, _, _ = multipletests(p_list, method='fdr_bh')
    df['adjusted_p_values'] = adjusted_p_values
    df.to_csv(output_file, index=False, sep='\t') 
"""

