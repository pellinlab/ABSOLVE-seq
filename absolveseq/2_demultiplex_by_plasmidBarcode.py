import pandas as pd
import numpy as np
import os
import re
import gzip
from tqdm import tqdm
import glob
import multiprocessing
import json
import ast
from collections import Counter
import argparse

def read_lines(file):
    header = file.readline().rstrip()
    seq = file.readline().rstrip()
    thirdline = file.readline().rstrip()
    qual = file.readline().rstrip()
    return header, seq, thirdline, qual

def write_lines(file, header, seq, thirdline, qual):
    file.write(header + os.linesep)
    file.write(seq + os.linesep)
    file.write(thirdline + os.linesep)
    file.write(qual + os.linesep)
    
def get_umis(r1):
    # may need to modify such that "ot" matches OT names in LVOT_OLIGO_POOL.xlsx
    dict = {} # save UMI, read, count info for sample
    r1_fa_fn = r1.split('/')[-1]
    # open FASTQ files for each OT to write to
    umi_list = []
    # read in FASTQ files for samples
    with gzip.open(r1, 'rt') as r1_fastq:
        while True:
            r1_header, r1_seq, r1_thirdline, r1_qual = read_lines(r1_fastq)
            umi_seq = r1_header.split(':')[-1]
            umi_list.append(umi_seq)
            if not r1_header: # EOF
                break
    return umi_list

def split_fastq_based_on_UMI_v2(r1_fn, amplicon_fn ='./test/data/target_info/OT_guide_amplicon_seq.csv', fa_out_folder='./singleUMI_fastq', 
                    crispress_input_folder='./singleUMI_crispresso_input', crispresso_output_folder="./CRISPResso_output", n_processes = 10):
    # if not os.path.exists(fa_out_folder):
    #     os.mkdir(fa_out_folder)
    # fa_out_folder='./singleUMI_fastq'
    os.makedirs(crispress_input_folder, exist_ok=True)
    os.makedirs(fa_out_folder, exist_ok=True)
    # crispress_input_folder = './singleUMI_crispresso_input'
    OT_name = r1_fn.split('/')[-1].split('_')[-2]
    sample_name = r1_fn.split('/')[-1].replace('_R1.fastq.gz', '')
    fa_out_folder = fa_out_folder + '/' + sample_name
    os.makedirs(fa_out_folder, exist_ok=True)
    all_umis = np.array(get_umis(r1_fn))[:-1]
    tdf = pd.DataFrame(index=all_umis)
    tdf['line_0'] = np.arange(len(all_umis))*4
    tdf = tdf.reset_index()
    tdf['line_0'] = tdf['line_0'].astype(str) 
    tdf = tdf.groupby('index')['line_0'].agg(','.join).reset_index()
    
    amplicon_df = pd.read_csv(amplicon_fn, index_col=['OT-Name'])
    amplicon_seq = amplicon_df.loc[OT_name, 'amplicon_seq']
    guide_seq = amplicon_df.loc[OT_name, 'guide_seq']
    cmd = ''
    input_df = pd.DataFrame()
    with gzip.open(r1_fn, 'rt') as r1_fastq,\
                gzip.open(r1_fn.replace("R1", "R2"), 'rt') as r2_fastq:
        r1_lines = r1_fastq.readlines()
        r2_lines = r2_fastq.readlines()
        for idx, row in tdf.iterrows():
            umi = row['index']
            select_line = row['line_0'].split(',')
            umi_r1_fn = fa_out_folder + '/' + umi + '_R1.fastq.gz' 
            # mode = 'at' if os.path.exists(umi_r1_fn) else 'wt'  # Append if exists, write if not
            ot_r1_lines = []
            ot_r2_lines = []
            for li in select_line:
                li = int(li)
                for i in range(4):
                    ot_r1_lines.append(r1_lines[li+i])
                    ot_r2_lines.append(r2_lines[li+i])
            ot_r1 = gzip.open(umi_r1_fn, 'wt')
            ot_r2 = gzip.open(umi_r1_fn.replace('R1', 'R2'), 'wt')
            ot_r1.writelines(ot_r1_lines)
            ot_r2.writelines(ot_r2_lines)
            ot_r1.close()
            ot_r2.close()
            # build input file for CRISPResso
            info = {'name':umi, 'fastq_r1':umi_r1_fn, 'fastq_r2':umi_r1_fn.replace('R1', 'R2'), 'amplicon_seq':amplicon_seq, 'guide_seq':guide_seq}
            input_df = pd.concat([input_df, pd.DataFrame(info, index=[0])])
            input_df = input_df[['name', 'fastq_r1', 'fastq_r2', 'amplicon_seq', 'guide_seq']]
            # if not os.path.exists(crispress_input_folder):
            #     os.mkdir(crispress_input_folder)
        input_fn = crispress_input_folder + '/' + sample_name + '_input.tsv'
        input_df.to_csv(input_fn,  index=False, sep='\t')
        # make crispressoBatch cmd
        cmd = "CRISPRessoBatch --batch_settings " + input_fn + " --batch_output_folder " + crispresso_output_folder + "/CRISPResso_on_" + sample_name + \
                    " --n_processes " +str(n_processes)+ " --ignore_substitutions --min_frequency_alleles_around_cut_to_plot 0 --skip_failed --write_detailed_allele_table --plot_window_size 10"
    return cmd

def main():
    parser = argparse.ArgumentParser(
        description='Split FASTQ files based on UMI in header and prepare CRISPResso input.'
    )
    parser.add_argument(
        '--fastq_dir',
        type=str,
        required=True,
        help='Directory containing FASTQ files with UMIs in the header.'
    )
    parser.add_argument(
        '--fa_out_folder',
        type=str,
        default='./test/demultiplexed_pBC_fastq',
        help='Output folder for single UMI FASTQ files.'
    )
    parser.add_argument(
        '--crispresso_input_folder',
        type=str,
        default='./test/CRISPResso_input_files',
        help='Output folder for CRISPResso input files.'
    )
    parser.add_argument(
        '--crispresso_output_folder',
        type=str,
        default='./CRISPResso_output',
        help='Output folder for CRISPResso results.'
    )
    parser.add_argument(
        '--amplicon_fn',
        type=str,
        default='./test/data/target_info/OT_guide_amplicon_seq.csv',
        help='CSV file containing amplicon and guide sequences.'
    )
    parser.add_argument(
        '--n_processes',
        type=int,
        default=10,
        help='Number of worker processes to use in parallel.'
    )
    
    args = parser.parse_args()

    r1_fn_list = glob.glob(os.path.join(args.fastq_dir, '*_R1.fastq.gz'))
    print('Total samples:', len(r1_fn_list))

    cmd_list = []

    # parallel over r1_fn_list
    pool_args = [
        (r1_fn, args.amplicon_fn, args.fa_out_folder, args.crispresso_input_folder, args.crispresso_output_folder, args.n_processes)
        for r1_fn in r1_fn_list
    ]

    # use multiprocessing.Pool for parallelism
    with multiprocessing.Pool(processes=args.n_processes) as pool:
        for cmd in tqdm(
            pool.starmap(split_fastq_based_on_UMI_v2, pool_args),
            total=len(pool_args),
            desc="Processing samples",
            unit="sample",
            dynamic_ncols=True
        ):
            cmd_list.append(cmd)

    with open('./CRISPRessoBatch_absolveseq.sh', 'w') as f:
        for cmd in cmd_list:
            f.write(cmd + '\n')

    print('CRISPRessoBatch commands saved to CRISPRessoBatch_absolveseq.sh')

if __name__ == "__main__":
    main()


