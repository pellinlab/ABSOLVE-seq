import pandas as pd
import gzip
import glob
import numpy as np
import os
from tqdm import tqdm
import argparse
import multiprocessing as mp  # <-- add this
# import seaborn as sns

def read_lines(file):
    header = file.readline().rstrip()
    seq = file.readline().rstrip()
    thirdline = file.readline().rstrip()
    qual = file.readline().rstrip()
    return header, seq, thirdline, qual

def write_lines(file, header, seq, thirdline, qual):
    with gzip.open(file, 'at') as file:
        file.write(header + os.linesep)
        file.write(seq + os.linesep)
        file.write(thirdline + os.linesep)
        file.write(qual + os.linesep)

def put_umi_in_header(r1_info, files_UMIs, OLIGO_INFO, n_thr=5):
    # r1 = "./wtCas9_trimmed/trimmed_Nola_EP031022_2_LVpool_wtCas9_D3734_T1_R1_paired.fastq.gz"
    # print(r1)
    total_read_count = 0
    filtered_read_count = 0 # reads that contain oligo barcode
    dict = {} # UMI counts
    OLIGO_POOL = OLIGO_INFO.set_index("Barcode")["OT"].to_dict()
    r1 = r1_info['Fastq_path']
    with gzip.open(r1, 'rt') as r1_fastq,\
        gzip.open(r1.replace("R1","R2"), 'rt') as r2_fastq:
        donor = r1_info['donor_id']
        rep = r1_info['replicate']
        cas_type = r1_info['cas_type']
        # sample = os.path.basename(r1)
        sample = '_'.join([cas_type, donor, rep])
        while True:
            r1_header, r1_seq, r1_thirdline, r1_qual = read_lines(r1_fastq)
            r2_header, r2_seq, r2_thirdline, r2_qual = read_lines(r2_fastq)
            if not r1_header: # EOF
                break
            oligo_barcode = r1_seq[39:48] # target ID
            constant_seq = r1_seq[16:16+23]
            total_read_count += 1
            # constant_seq 
            pass_flag = True
            n_count = 0
            expected_const =  "CCAACCTCATAGAACACTCATCC"
            for i, nt in enumerate(constant_seq):
                if nt == "N":
                    n_count += 1
                elif nt != expected_const[i]:
                    pass_flag = False
                    break
            if n_count > n_thr:
                pass_flag = False
            # if valid barcode and correct UMI length
            if (oligo_barcode in OLIGO_POOL) and (pass_flag):
                ot = OLIGO_POOL[oligo_barcode]
                ot = '_'.join([cas_type, donor, rep, ot])
                # get FASTQs to write to 
                ot_r1 = files_UMIs[ot][0]
                ot_r2 = files_UMIs[ot][1]
                # extract UMI from read 1
                umi_seq = r1_seq[:16]
                # add UMI to header (immediately after last ":")
                r1_header_umi = r1_header + ":" + umi_seq
                r2_header_umi = r2_header + ":" + umi_seq
                # remove UMI from both read 1 & read 2
                r1_no_umi_seq = r1_seq[16:] # remove first 16 nt (UMI) from read 1
                r1_qual = r1_qual[16:]
                r2_no_umi_seq = r2_seq[:-16] # remove last 16 nt (UMI) from read 2
                r2_qual = r2_qual[:-16]
                # write to new FASTQ files based on OT corresponding to oligo barcode
                write_lines(ot_r1, r1_header_umi, r1_no_umi_seq, r1_thirdline, r1_qual)
                write_lines(ot_r2, r2_header_umi, r2_no_umi_seq, r2_thirdline, r2_qual)
                filtered_read_count += 1
                # save UMI & count info for sample
                if umi_seq in dict:
                    dict[umi_seq] += 1
                else:
                    dict[umi_seq] = 1
                # save UMIs per OT for samples
                if umi_seq in files_UMIs[ot][2]:
                    files_UMIs[ot][2][umi_seq] += 1
                else:
                    files_UMIs[ot][2][umi_seq] = 1
    # close files and count UMIs per OT
    UMIs_per_OT = [] # for sample
    for ot in files_UMIs:
        UMIs_per_OT.append(len(files_UMIs[ot][2]))

    # return saved info for sample
    info = [sample, total_read_count, filtered_read_count, len(dict)]
    return dict, info, UMIs_per_OT

def _process_one_row(args):
    """
    Worker wrapper so it can be pickled by multiprocessing.
    Builds files_UMIs and calls put_umi_in_header for one sample row.
    """
    row_dict, out_dir, oligo_info_dict = args
    row = pd.Series(row_dict)
    OLIGO_INFO = pd.DataFrame(oligo_info_dict)  # reconstruct DataFrame in worker
    OLIGO_POOL = OLIGO_INFO.set_index("Barcode")["OT"].to_dict()

    donor = row['donor_id']
    rep = row['replicate']
    cas_type = row['cas_type']
    sample = '_'.join([cas_type, donor, rep])

    files_UMIs = {}
    for ot in OLIGO_POOL.values():
        ot_r1 = os.path.join(out_dir, f"{sample}_{ot}_R1.fastq.gz")
        ot_r2 = os.path.join(out_dir, f"{sample}_{ot}_R2.fastq.gz")
        files_UMIs[f"{sample}_{ot}"] = [ot_r1, ot_r2, {}]
        # touch files so write_lines() can open in 'at'
        with gzip.open(ot_r1, 'wt'):
            pass
        with gzip.open(ot_r2, 'wt'):
            pass

    return put_umi_in_header(row, files_UMIs, OLIGO_INFO)

def main():
    parser = argparse.ArgumentParser(description="Extract plasmid barcodes and move UMIs to read headers.")
    parser.add_argument(
        '--output_dir',
        type=str,
        default='./test/demultiplexed_tBC_fastq',
        help='Directory to save output FASTQ files and summaries.'
    )
    parser.add_argument(
        '--target_info',
        type=str,
        default='./test/data/target_info/LVOT_oligo_pool.xlsx',
        help='Excel file containing oligo pool information.'
    )
    parser.add_argument(
        '--sample_info',
        type=str,
        default='./test/data/target_info/NovaSeq3_sample_info_example.csv',
        help='csv file containing oligo pool information.'
    )
    parser.add_argument(
        '--n_cpu',
        type=int,
        default=4,
        help='Number of worker processes.'
    )

    args = parser.parse_args()
    out_dir = args.output_dir
    target_info = args.target_info
    sample_info = args.sample_info
    n_cpu = args.n_cpu

    os.makedirs(out_dir, exist_ok=True)
    OLIGO_INFO = pd.read_excel(target_info, engine = "openpyxl")
    OLIGO_POOL = OLIGO_INFO.set_index("Barcode")["OT"].to_dict()

    sample_df = pd.read_csv(sample_info)

    # prepare arguments for workers
    oligo_info_dict = OLIGO_INFO.to_dict(orient="list")
    worker_args = [
        (row.to_dict(), out_dir, oligo_info_dict)
        for _, row in sample_df.iterrows()
    ]

    # parallel per-sample processing
    with mp.Pool(processes=n_cpu) as pool:
        results = list(
            tqdm(pool.imap_unordered(_process_one_row, worker_args),
                 total=len(worker_args))
        )

    dict_list = []
    info_list = []
    UMIs_per_OT_list = []
    for d, info, umis in results:
        dict_list.append(d)
        info_list.append(info)
        UMIs_per_OT_list.append(umis)

    per_ot_df = pd.DataFrame.from_dict(OLIGO_POOL, orient='index')
    per_ot_df.reset_index(drop=True, inplace=True)

    # for each sample
    for i in range(len(info_list)):
        sample = info_list[i][0]
        # transform dict into df
        df = pd.DataFrame(dict_list[i], index=[0]).T

        # save file with detailed UMI info
        df = df.sort_values(by=0, ascending=False)
        df.index.rename("plasmid_barcode", inplace=True)
        df.columns = ["Reads"]
        os.makedirs(os.path.join(out_dir, "pBCs_stats"), exist_ok=True)
        df.to_csv(os.path.join(out_dir, "pBCs_stats", f"{sample}.csv"))
        # add per OT UMI counts
        per_ot_df[sample] = UMIs_per_OT_list[i]

    # write read & UMI counts summaries to file
    summary_df = pd.DataFrame(
        info_list,
        columns=["Sample","Total reads", "Reads with plasmid barcode", "#plasmid barcodes"]
    )
    with pd.ExcelWriter(os.path.join(out_dir, "absovleseq_reads_pBCs.xlsx")) as writer:
        summary_df.to_excel(writer, sheet_name="Per_sample_summary", index=False)
        per_ot_df.to_excel(writer, sheet_name="Per_OT_summary", index=False)

if __name__ == "__main__":
    main()