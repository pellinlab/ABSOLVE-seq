import pandas as pd
import numpy as np
import os
import re
import gzip
from tqdm import tqdm
from multiprocessing import Pool
import glob
import multiprocessing
import json
import ast
import argparse
import math
from collections import Counter

# Functions to process UMI alleles
def filter_low_reads_umiallele(df):
    filtered_df = df[df['#Reads'] > 1].reset_index()
    filtered_df['%Reads'] = filtered_df['#Reads']/filtered_df['#Reads'].sum() * 100
    return filtered_df

def deduplicate_umiallele(df, filtering=False):
    tdf = df.sort_values(by=['UMI', '%Reads_UMI'], ascending=False)
    tdf_dedup = tdf.drop_duplicates(subset='UMI', keep='first')
    if filtering:
        filtered_df = tdf_dedup[tdf_dedup['%Reads_UMI'] > 65].reset_index()
    else:
        filtered_df = tdf_dedup.reset_index()
    filtered_df['#Reads'] = 1
    filtered_df['%Reads'] = filtered_df['#Reads']/filtered_df['#Reads'].sum() * 100
    return filtered_df

def levenshteinFullMatrix(str1, str2):
    m = len(str1)
    n = len(str2)
    # Initialize a matrix to store the edit distances
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize the first row and column with values from 0 to m and 0 to n respectively
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    # Fill the matrix using dynamic programming to compute edit distances
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                # Characters match, no operation needed
                dp[i][j] = dp[i - 1][j - 1]
            else:
                # Characters don't match, choose minimum cost among insertion, deletion, or substitution
                dp[i][j] = 1 + min(dp[i][j - 1], dp[i - 1][j], dp[i - 1][j - 1])

    # Return the edit distance between the strings
    return dp[m][n]

def barcode_entropy(sequence):
    """
    Calculate the Shannon entropy of a DNA sequence.
    """
    # Count the frequency of each nucleotide
    counts = Counter(sequence)
    total = len(sequence)
    
    # Calculate the entropy
    entropy_value = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy_value

#  add parameter
argparser = argparse.ArgumentParser()
argparser.add_argument('--plasmid_barcode_annot_file', type=str, help='Path to the plamsid barcode annotation file', default='./test/data/plasmid_barcode_category.tsv.gz')
argparser.add_argument('--crispresso_alleles_dir',  type=str, help='Path to the preprocesed editing allele tables from crispresso2', default='./test/absolveseq_edits/crispresso_allele_tables/')
argparser.add_argument('--output_dir',  type=str, help='annotated allele table directory', default='./test/absolveseq_edits/raw/')
maxi_ds_df = pd.read_csv(argparser.parse_args().plasmid_barcode_annot_file, sep='\t')
editing_alleles_dir = argparser.parse_args().crispresso_alleles_dir
out_dir = argparser.parse_args().output_dir
os.makedirs(out_dir, exist_ok=True)
# print(maxi_ds_df.columns)
# allele_table_dir = './alleles_HiFi_WTCas9_NoEP/'
for fn in tqdm(glob.glob(editing_alleles_dir + '/*Allele_frequency.txt')):
    ot_name = fn.split('/')[-1].split('_')[0].replace('-', '_')
    # if ot_name not in sel_OTs:
    #     continue
    out_name = fn.split('/')[-1]
    allele_df = pd.read_csv(fn, sep='\t')
    # print(allele_df.columns)
    # ot_name = allele_df.iloc[0]['OT_name'].replace('-', '_')
    OT_maxi_df = maxi_ds_df[(maxi_ds_df['target_matched'] == ot_name)].drop_duplicates(subset='plasmid_barcode', keep='first')
    # umi_set = list(set(allele_df['UMI'].unique())
    # allele_df = allele_df.merge(umi_complexity.loc[umi_set][['#barcode', 'entropy']], left_on='UMI', right_index=True, how='left')
    allele_maxi_df = allele_df.merge(OT_maxi_df, how='left', left_on='UMI', right_on='plasmid_barcode')
    # allele_maxi_df.loc[allele_maxi_df[allele_maxi_df['good_L0'].isna()].index, 'good_L0'] = 'Unseen'
    # allele_maxi_df.loc[allele_maxi_df[allele_maxi_df['with_expected_OT'].isna()].index, 'with_expected_OT'] = 'Unseen'
    # allele_maxi_df['good_L0'] = allele_maxi_df['good_L0'].astype(str)
    # Stringent
    allele_maxi_df.loc[allele_maxi_df[allele_maxi_df['if_good'].isna()].index, 'if_good'] = 'Unseen'
    # allele_maxi_df.loc[allele_maxi_df[allele_maxi_df['with_expected_OT'].isna()].index, 'with_expected_OT'] = 'Unseen'
    allele_maxi_df['if_good'] = allele_maxi_df['if_good'].astype(str)
    # print(allele_maxi_df.columns)
    # allele_maxi_df['with_expected_OT'] = allele_maxi_df['with_expected_OT'].astype(str)
    allele_maxi_df.to_csv(out_dir + '/' + out_name.replace('.txt', '.annot.txt'), index=False, sep='\t')
    break