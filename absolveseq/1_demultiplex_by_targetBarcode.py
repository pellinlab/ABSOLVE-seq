"""
Linda Lin
4/13/2022

Input: LV OT PE FASTQs post trimming (with UMIs within reads) for pool

Output:
- LV OT PE FASTQs demultiplexed by oligo barcode (filter out reads without expected barcode and/or unexpected UMI length), with UMIs moved to read1/read2 header (at end after last ":")
- Excel spreadsheet with summary of read & UMI counts by sample and by OT
- CSV files for each sample with UMIs

Notes:
- May need to rename samples and/or modify code (lines 28-30) rename samples such that "ot" matches OT names in LVOT_OLIGO_POOL.xlsx
- Run script in folder with FASTQs present as well as empty dirs "umi_in_header/umis"
- Can process many samples in parallels
"""

import gzip
import os
import glob
import multiprocessing as mp
import pandas as pd
import argsparse



def put_umi_in_header(r1, constant_seq="CCAACCTCATAGAACA"):
    # may need to modify such that "ot" matches OT names in LVOT_OLIGO_POOL.xlsx
    sample = r1.split("_")[1]
    ot = sample[21:]

    total_read_count = 0
    filtered_read_count = 0 # reads that contain oligo barcode
    dict = {} # save UMI, read, count info for sample
    files_UMIs = {} # save file names and UMIs per OT

    # open FASTQ files for each OT to write to
    for ot in OLIGO_POOL.values():
        ot_r1 = gzip.open(out_dir + "/" + sample + "_" + ot + "_R1.fastq.gz", 'wt')
        ot_r2 = gzip.open(out_dir + "/" + sample + "_" + ot + "_R2.fastq.gz", 'wt')
        files_UMIs[ot] = [ot_r1, ot_r2, {}]

    # read in FASTQ files for samples
    with gzip.open(r1, 'rt') as r1_fastq,\
         gzip.open(r1.replace("R1","R2"), 'rt') as r2_fastq:

        while True:
            r1_header, r1_seq, r1_thirdline, r1_qual = read_lines(r1_fastq)
            r2_header, r2_seq, r2_thirdline, r2_qual = read_lines(r2_fastq)
            if not r1_header: # EOF
                break

            # extract oligo barcode
            oligo_barcode = r1_seq[39:48]
            total_read_count += 1

            # get start of constant sequence
            umi_end = r1_seq.find(constant_seq)

            # if valid barcode and correct UMI length
            if (oligo_barcode in OLIGO_POOL) and (umi_end == 16):
                ot = OLIGO_POOL[oligo_barcode]
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
        files_UMIs[ot][0].close()
        files_UMIs[ot][1].close()
        UMIs_per_OT.append(len(files_UMIs[ot][2]))

    # return saved info for sample
    info = [sample, total_read_count, filtered_read_count, len(dict)]
    return dict, info, UMIs_per_OT


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


parser = argsparse.ArgumentParser(description="Extract plasmid barcodes and move UMIs to read headers.")
parser.add_argument('--fastq_dir', type=str, default='.', help='Directory containing input FASTQ files.')
parser.add_argument('--output_dir', type=str, default='../test/fastq_with_pBC', help='Directory to save output FASTQ files and summaries.')

out_dir = parser.parse_args().output_dir
# expected oligo barcode - OT correspondences
OLIGO_POOL = pd.read_excel("../test/data/target_info/LVOT_oligo_pool.xlsxLVOT_oligo_pool.xlsx", engine = "openpyxl")
OLIGO_POOL = OLIGO_POOL.set_index("Barcode")["OT"].to_dict()

pool = mp.Pool(10)
fastq_dir = parser.parse_args().fastq_dir
fastq_files = glob.glob(os.path.join(fastq_dir, "*_R1_*"))

# process each sample
dict, info, UMIs_per_OT = zip(*pool.map(put_umi_in_header, fastq_files))

pool.close()
pool.join()

umi_df = pd.DataFrame()
per_ot_df = pd.DataFrame.from_dict(OLIGO_POOL, orient='index')
per_ot_df.reset_index(drop=True)
os.makedirs(out_dir + "/umis/", exist_ok=True)
# for each sample
for i in range(len(info)):
    sample = info[i][0]

    # transform dict into df
    df = pd.DataFrame(dict[i], index=[0]).T

    # save file with detailed UMI info
    df = df.sort_values(by=0, ascending=False)
    df.index.rename("UMI", inplace=True)
    df.columns = ["Reads"]
    df.to_csv(out_dir + "/umis/" + sample + ".csv")
    # add per OT UMI countsa
    per_ot_df[sample] = UMIs_per_OT[i]

# write read & UMI counts summaries to file
summary_df = pd.DataFrame(info, columns=["Sample","Total reads", "Reads with oligo barcode", "Corresponding UMIs"])
with pd.ExcelWriter(out_dir + "/LVOT_reads_UMIs.xlsx") as writer:
    summary_df.to_excel(writer, sheet_name="Per_sample_summary", index=False)
    per_ot_df.to_excel(writer, sheet_name="Per_OT_summary", index=False)