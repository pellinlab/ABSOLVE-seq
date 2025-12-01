import glob
import os
import pandas as pd
from tqdm import tqdm
import json
import argparse

def add_cutsite(allele_full_df):
    # out_name = fn.split('/')[-1]
    # allele_full_df = pd.read_csv(fn, sep='\t')
    align_20bp = []
    ref_20bp = []
    for idx, row in allele_full_df.iterrows():
        ref_pos = np.array(ast.literal_eval(row['ref_positions']))
        cutsite = row['sgRNA_cut_points']
        cutsite = np.arange(len(row['Reference_Sequence']))[ref_pos == cutsite][0]
        align_20bp.append(row['Aligned_Sequence'][cutsite-9:cutsite+11])
        ref_20bp.append(row['Reference_Sequence'][cutsite-9:cutsite+11])
    allele_full_df['Unedited'] = allele_full_df['Read_Status'] == 'UNMODIFIED'
    allele_full_df['Aligned_Sequence_20bp'] = align_20bp
    allele_full_df['Reference_Sequence_20bp'] = ref_20bp
    allele_full_df.to_csv(fn, sep='\t', index=False)
    return True

def combine_allele_table(sample_name, crisprsso_res_folder = './test/CRISPResso_output/', out_folder = './absolveseq_edits/alleles_all/'):
    os.makedirs(out_folder, exist_ok=True)
    ot = sample_name.split('_')[-1]
    group_i = sample_name.split('_')[0]
    donor_i = sample_name.split('_')[1]
    Rep_i = sample_name.split('_')[2]

    out_name = ot + '_' + donor_i + '_' + group_i + '_' + Rep_i
    alleles_dir = crisprsso_res_folder +'/CRISPResso_on_'+sample_name+'/CRISPRessoBatch_on_'+sample_name+'_input/CRISPResso_on_*/'
    print(f"Alleles directory: {alleles_dir}")
    alleles_fns = glob.glob(crisprsso_res_folder +'/CRISPResso_on_'+sample_name+'/CRISPRessoBatch_on_'+sample_name+'_input/CRISPResso_on_*/Alleles_frequency_table.zip')
    umi_full_allele_list = []
    for fn in alleles_fns:
        umi_seq = fn.split('/')[-2].split('_')[-1]
        full_df = pd.read_csv(fn, sep='\t')
        json_file = fn.replace(fn.split('/')[-1], 'CRISPResso2_info.json')
        json_data = json.load(open(json_file))
        cutsite = json_data['results']['refs']['Reference'] ['sgRNA_cut_points'][0]
        full_df['UMI'] = umi_seq
        full_df['sgRNA_cut_points'] = cutsite
        umi_full_allele_list.append(full_df)
    umi_full_allele_df = pd.concat(umi_full_allele_list).reset_index()
    umi_full_allele_df = umi_full_allele_df.rename(columns={'%Reads':'%Reads_UMI'})
    umi_full_allele_df['%Reads'] = umi_full_allele_df['#Reads']/umi_full_allele_df['#Reads'].sum() * 100
    umi_full_allele_df['sample_id'] = sample_name
    umi_full_allele_df['donor'] = donor_i
    umi_full_allele_df['replicate'] = Rep_i
    umi_full_allele_df['group'] = group_i
    umi_full_allele_df['OT_name'] = ot
    umi_full_allele_df = add_cutsite(umi_full_allele_df)
    umi_full_allele_df = umi_full_allele_df[['UMI', 'index', 'Aligned_Sequence',
               'Reference_Sequence', 'Aligned_Sequence_20bp', 'Reference_Sequence_20bp', 'sgRNA_cut_points', 'Aligned_Reference_Scores', 'Unedited', 'n_deleted', 'n_inserted',
               'n_mutated', '#Reads', '%Reads_UMI', '%Reads', 'sample_id', 'donor',
               'replicate', 'group', 'OT_name']]
    umi_full_allele_df.to_csv(out_folder + '/' + out_name + '_Allele_frequency_table_full_withUMI.txt', index=False, sep='\t')
    return out_name

argparser = argparse.ArgumentParser()
argparser.add_argument('--crispresso_result_dir', type=str, help='Path to the CRISPResso results folder', default='./test/CRISPResso_output/')
argparser.add_argument('--out_folder', type=str, help='Path to the output folder', default='./absolveseq_edits/crispresso_allele_tables/')
crispresso_out_folder = argparser.parse_args().crispresso_result_dir
out_folder = argparser.parse_args().out_folder
os.makedirs(out_folder, exist_ok=True)

# collect alleles table and clean the /dev/shm outputs
for ot_fn in tqdm(glob.glob(crispresso_out_folder + '/*')):
    basename = os.path.basename(ot_fn).replace('CRISPResso_on_', '')
    combine_allele_table(basename, out_folder = out_folder, crisprsso_res_folder=crispresso_out_folder)
    # print('combine_allele_table finished')
    # break