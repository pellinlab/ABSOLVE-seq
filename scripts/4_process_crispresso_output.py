import glob
import os
import pandas as pd
import json

def combine_allele_table(sample_name, out_folder = '/cpl/home/linjc/project/ABSOLVE-seq/singleUMI_crispresso_WTCas9/alleles_all/', return_df = False):
    os.makedirs(out_folder, exist_ok=True)
    # lvid = sample_name.split('_')[0].replace('-', '_')
    LV_cvt = pd.read_csv('/cpl/home/linjc/project/ABSOLVE-seq/metadata/NovaSeq3_sample_info_new.csv').set_index('name')
    # donor_i + '_' + group_i + '_' + Rep_i
    # LV_cvt['name'] = LV_cvt['donor_id'] + '_' + LV_cvt['cas_type'] + '_' + LV_cvt['replicate']
    # LV_cvt = LV_cvt.set_index('name')
    ot = sample_name.split('_')[-1]
    donor_i = sample_name.split('_')[1]
    group_i = sample_name.split('_')[0]+'Cas9'
    Rep_i = sample_name.split('_')[2]
    name = donor_i + '_' + group_i + '_' + Rep_i
    lvid = LV_cvt.loc[name, 'LV_ID']
    # print(lvid, ot)
    # LV_cvt = pd.read_csv('./metadata/NovaSeq3_sample_info.csv')
    # Rep_i =  LV_cvt.loc[lvid, 'replicate']
    # group_i = LV_cvt.loc[lvid, 'cas_type'].replace('Cas9', '').replace('noEP', 'NoEP')
    # donor_i = LV_cvt.loc[lvid, 'donor_id']
    out_name = ot + '_' + lvid.replace('_', '-') + '_'+ donor_i + '_' + group_i + '_' + Rep_i
    # if not os.path.exists(out_folder):
    #     os.mkdir(out_folder)
    # CRISPResso_on_LV-S1_1450-OT-0000-REF/CRISPRessoBatch_on_LV-S1_1450-OT-0000-REF_input
    alleles_fns = glob.glob('/dev/shm/jiecong/singleUMI_crispresso/CRISPResso_on_'+sample_name+'/CRISPRessoBatch_on_'+sample_name+'_input/CRISPResso_on_*/Alleles_frequency_table_around_sgRNA*')
    # print(len(alleles_fns))
    umi_allele_list = []
    umi_full_allele_list = []
    for fn in alleles_fns:
        umi_seq = fn.split('/')[-2].split('_')[-1]
        df = pd.read_csv(fn, sep='\t')
        full_df = pd.read_csv(fn.replace(fn.split('/')[-1], 'Alleles_frequency_table.zip'), sep='\t')
        json_file = fn.replace(fn.split('/')[-1], 'CRISPResso2_info.json')
        json_data = json.load(open(json_file))
        cutsite = json_data['results']['refs']['Reference'] ['sgRNA_cut_points'][0]
        # json_data['']
        # new_json_fn = fn.replace('CRISPResso2_info.json', '{}_CRISPResso2_info.json',)
        # os.move(json_file, json_file.replace('CRISPResso2_info.json', '{}_CRISPResso2_info.json'.format(umi_seq)))
        full_df['UMI'] = umi_seq
        full_df['sgRNA_cut_points'] = cutsite
        df['UMI'] = umi_seq
        umi_allele_list.append(df)
        umi_full_allele_list.append(full_df)
    umi_allele_df = pd.concat(umi_allele_list).reset_index()
    umi_full_allele_df = pd.concat(umi_full_allele_list).reset_index()
    umi_allele_df = umi_allele_df.rename(columns={'%Reads':'%Reads_UMI'})
    umi_full_allele_df = umi_full_allele_df.rename(columns={'%Reads':'%Reads_UMI'})
    umi_allele_df['%Reads'] = umi_allele_df['#Reads']/umi_allele_df['#Reads'].sum() * 100
    umi_full_allele_df['%Reads'] = umi_full_allele_df['#Reads']/umi_full_allele_df['#Reads'].sum() * 100
    # umi_allele_df['UMI|Alelle'] = umi_allele_df['UMI'] + '|' + umi_allele_df['Aligned_Sequence']
    umi_allele_df['sample_id'] = sample_name
    umi_full_allele_df['sample_id'] = sample_name
    umi_allele_df['donor'] = donor_i
    umi_full_allele_df['donor'] = donor_i
    umi_allele_df['replicate'] = Rep_i
    umi_full_allele_df['replicate'] = Rep_i
    umi_allele_df['group'] = group_i
    umi_full_allele_df['group'] = group_i
    umi_allele_df['OT_name'] = ot
    umi_full_allele_df['OT_name'] = ot
    umi_allele_df['LV_id'] = lvid.replace('_', '-')
    umi_full_allele_df['LV_id'] = lvid.replace('_', '-')
    # print(out_name)
    if return_df:
        return umi_allele_df, umi_full_allele_df
    else:
        umi_allele_df.to_csv(out_folder + '/' + out_name + '_Allele_frequency_table_withUMI.txt', index=False, sep='\t')
        umi_full_allele_df.to_csv(out_folder + '/' + out_name + '_Allele_frequency_table_full_withUMI.txt', index=False, sep='\t')
        return out_name

# collect alleles table and clean the /dev/shm outputs
files = ['MODIFICATION_FREQUENCY_SUMMARY.txt', 'MODIFICATION_PERCENTAGE_SUMMARY.txt', 'Nucleotide_frequency_summary.txt', 'Nucleotide_percentage_summary.txt']
crispresso_out_folder = '/dev/shm/jiecong/singleUMI_crispresso/'
for ot_fn in glob.glob(crispresso_out_folder + '*'):
    basename = os.path.basename(ot_fn).replace('CRISPResso_on_', '')
    # print(basename)
    input_folder = os.path.join(ot_fn, 'CRISPRessoBatch_on_' + basename + '_input/') 
    print(input_folder)
    finished_flag = True
    for finished_fn in files:
        finished_fn = os.path.join(input_folder, finished_fn)
        if not os.path.exists(finished_fn):
            finished_flag = False
            print(basename, 'not finished')
            break
    if finished_flag:
        print(basename, 'finished')
        combine_allele_table(basename, out_folder = '/cpl/home/linjc/project/ABSOLVE-seq/singleUMI_crispresso_WTCas9/alleles_all/')
        print('combine_allele_table finished')
        # os.system('tar -czf {}.tar.gz {} > /dev/null 2>&1'.format(input_folder, input_folder))
        # os.system('rm -rf {}'.format(ot_fn))
        print(ot_fn, 'deleted orignial folder')