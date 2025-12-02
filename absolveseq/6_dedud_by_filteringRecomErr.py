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


def deduplicate_umiallele(df, filtering=False):
    tdf = df.sort_values(by=['UMI', '%Reads_UMI'], ascending=False)
    tdf_dedup = tdf.drop_duplicates(subset='UMI', keep='first')
    if filtering:
        filtered_df = tdf_dedup[tdf_dedup['%Reads_UMI'] > 65].reset_index()
    else:
        filtered_df = tdf_dedup.reset_index()
    filtered_df['#Reads'] = 1
    filtered_df['%Reads'] = filtered_df['#Reads'] / filtered_df['#Reads'].sum() * 100
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

# argparse.ArgumentParser()
parser = argparse.ArgumentParser()
parser.add_argument('--target_oligo_file', type=str, help='Path to the target oligo sequences file', default='./test/data/target_info/target_oligos_sequenes.csv')
parser.add_argument('--dedud_alleles_dir', type=str, help='Path to the processed editing allele tables from CRISPResso2', default='./test/absolveseq_edits/dedud/')
parser.add_argument('--output_dir', type=str, help='Directory for filtered deduplicated allele tables', default='./test/absolveseq_edits/dedud_filtered/')
args = parser.parse_args()
oligo_df = pd.read_csv(args.target_oligo_file).set_index('OT')
dir_path = args.dedud_alleles_dir
out_folder = args.output_dir
# for filter in ['dedud']']:
for fn in tqdm(glob.glob(dir_path + '/*.annot.txt.gz')):
    out_name = fn.split('/')[-1]
    ot_name = out_name.split('_')[0]
    barcode = oligo_df.loc[ot_name, 'Barcode']
    constant = 'CCAACCTCATAGAACACTCATCC'
    target_seq_pos = len(barcode) + len(constant)
    allele_df = pd.read_csv(fn, sep='\t')
    # allele_df = allele_df.drop(columns=['best_match_strigent'])
    allele_df['best_match'] = ot_name
    allele_df['same_match'] = ''
    allele_df['#same_match'] = 0
    allele_edited = allele_df[allele_df['Unedited'] == False]
    recomb_edits_idx = []
    mis_targets = []
    mis_targets_score = []
    same_edits_idx = []
    same_targets = [] 
    same_targets_score = []
    n_same_targets = []
    for idx, row in allele_edited.iterrows():
        target_seq = row['Aligned_Sequence'].replace('-', '')
        target_seq = target_seq[target_seq_pos:]
        expected_dist = levenshteinFullMatrix(target_seq, oligo_df.loc[ot_name, 'OT_sequence'])
        current_dist = expected_dist
        # print(ot_name, target_seq, oligo_df.loc[ot_name.replace('-', '_'), 'OT_sequence'], current_dist)
        best_match = ot_name
        same_match = []
        # best_match_S = ot_name
        for ot, row in oligo_df.iterrows():
            ot = ot.replace('_', '-')
            if ot == ot_name:
                continue
            tdist = levenshteinFullMatrix(target_seq, row['OT_sequence'])
            if tdist < current_dist:
                # print('\t', ot, tdist, target_seq, row['OT_sequence'])
                best_match = ot
                current_dist = tdist
            if tdist == expected_dist:
                same_match.append(ot)

        if best_match != ot_name:
            recomb_edits_idx.append(idx)
            mis_targets.append(best_match)
            mis_targets_score.append(current_dist)
        if len(same_match) > 0:
            same_edits_idx.append(idx)
            same_targets.append(','.join(same_match))
            n_same_targets.append(len(same_match))
            same_targets_score.append(expected_dist)
            # print(ot_name, best_match)
    allele_df.loc[recomb_edits_idx, 'best_match'] = mis_targets
    allele_df.loc[same_edits_idx, 'same_match'] = same_targets
    allele_df.loc[same_edits_idx, '#same_match'] = n_same_targets
    allele_df.loc[recomb_edits_idx, 'best_match_score'] = mis_targets_score
    allele_df.loc[same_edits_idx, 'same_match_score'] = same_targets_score
    out_fn = out_folder + out_name
    allele_df.to_csv(out_fn, sep='\t', index=False, compression='gzip')

parent_dir = os.path.dirname(os.path.dirname(out_folder))
for fn in tqdm(glob.glob(dir_path + '/*.annot.txt.gz')):
    out_name = fn.split('/')[-1]
    allele_df = pd.read_csv(fn, sep='\t')
    allele_df = allele_df[(allele_df['OT_name'] == allele_df['best_match'])&(allele_df['#same_match'] == 0)].reset_index(drop=True)
    # if not os.path.exists(out_dir):
    #     os.mkdir(out_dir)
    # out_fn = out_dir + out_name
    allele_df.to_csv(fn, sep='\t', index=False)
    # allele_df['UMI-Aligned_Sequence'] = allele_df['UMI'] + '-' + allele_df['Aligned_Sequence']
    allele_deduplicate_df = deduplicate_umiallele(allele_df)# .drop(columns='level_0')
    out_folder = parent_dir + '/dedud_filtered_dedup/'
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    out_fn = out_folder + '/' + out_name
    allele_deduplicate_df.to_csv(out_fn, sep='\t', index=False)