import pandas as pd
import numpy as np
import os
import re
import gzip
from tqdm import tqdm
import glob
import argparse
import math
from collections import Counter

def filter_low_reads_umiallele(df):
    filtered_df = df[df['#Reads'] > 1].reset_index()
    filtered_df['%Reads'] = filtered_df['#Reads'] / filtered_df['#Reads'].sum() * 100
    return filtered_df

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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plasmid_barcode_annot_file', type=str, help='Path to the plasmid barcode annotation file', default='./test/data/plasmid_barcode_category.tsv.gz')
    parser.add_argument('--crispresso_alleles_dir', type=str, help='Path to the processed editing allele tables from CRISPResso2', default='./test/absolveseq_edits/crispresso_allele_tables/')
    parser.add_argument('--output_dir', type=str, help='Directory for annotated allele tables', default='./test/absolveseq_edits/')
    args = parser.parse_args()

    maxi_ds_df = pd.read_csv(args.plasmid_barcode_annot_file, sep='\t')
    editing_alleles_dir = args.crispresso_alleles_dir
    # out_dir = args.output_dir
    dedud_dir = args.output_dir
    parent_dir = os.path.dirname(os.path.dirname(dedud_dir))
    
    raw_out_dir = parent_dir + '/raw/'
    os.makedirs(raw_out_dir, exist_ok=True)
    # print(raw_out_dir)
    print("Annotating ABSOLVE-seq editing outcome allele tables at {}".format(raw_out_dir))
    for fn in tqdm(glob.glob(os.path.join(editing_alleles_dir, '*Allele_frequency.txt'))):
        ot_name = os.path.basename(fn).split('_')[0].replace('-', '_')
        out_name = os.path.basename(fn)
        allele_df = pd.read_csv(fn, sep='\t')
        OT_maxi_df = maxi_ds_df[maxi_ds_df['target_matched'] == ot_name].drop_duplicates(subset='plasmid_barcode', keep='first')
        allele_maxi_df = allele_df.merge(OT_maxi_df, how='left', left_on='UMI', right_on='plasmid_barcode')
        allele_maxi_df.loc[allele_maxi_df['if_good'].isna(), 'if_good'] = 'Unseen'
        allele_maxi_df['if_good'] = allele_maxi_df['if_good'].astype(str)
        out_path = os.path.join(raw_out_dir, out_name.replace('.txt', '.annot.txt.gz'))
        allele_maxi_df.to_csv(out_path, index=False, sep='\t', compression='gzip')

    print("Deduding ABSOLVE-seq editing outcome...")
    for fn in tqdm(glob.glob(raw_out_dir +'/*.annot.txt.gz')):
        out_name = os.path.basename(fn)
        allele_df = pd.read_csv(fn, sep='\t')
        allele_deduplicate_df = deduplicate_umiallele(allele_df)
        out_folder = os.path.join(parent_dir, 'raw_dedup')
        os.makedirs(out_folder, exist_ok=True)
        out_fn = os.path.join(out_folder, out_name)
        allele_deduplicate_df.to_csv(out_fn, sep='\t', index=False, compression='gzip')

        dedud_df = allele_df[allele_df['if_good'] == 'True'].reset_index()
        dedud_df['%Reads'] = dedud_df['#Reads'] / dedud_df['#Reads'].sum() * 100
        if 'level_0' in dedud_df.columns:
            dedud_df = dedud_df.drop(columns='level_0')
        dedud_filtered_df = dedud_df[(dedud_df['entropy'] > 1.477) & ((dedud_df['#barcode'] == 1) | (dedud_df['#barcode'] == 2))]
        out_folder = os.path.join(dedud_dir)
        os.makedirs(out_folder, exist_ok=True)
        dedud_filtered_df.to_csv(os.path.join(out_folder, out_name), sep='\t', index=False, compression='gzip')

        allele_deduplicate_df = deduplicate_umiallele(dedud_filtered_df)
        out_folder = os.path.join(parent_dir, 'dedud_dedup')
        os.makedirs(out_folder, exist_ok=True)
        out_fn = os.path.join(out_folder, out_name)
        allele_deduplicate_df.to_csv(out_fn, sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    main()