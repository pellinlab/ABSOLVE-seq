# absolveseq: The ABSOLVE-Seq Analysis Package


The absolveseq package implements our data preprocessing and analysis pipeline for ABSOLVE-Seq data.

### References

##### The original paper describing the ABSOLVE-Seq method:

Jiecong Lin, My Anh Nguyen, Linda Y. Lin, Jing Zeng, Archana Verma, Nola R. Neri, Lucas Ferreira da Silva, Adele Mucci, Scot Wolfe, Kit L Shaw, Kendell Clement, Christian Brendel, Luca Pinello, Danilo Pellin, Daniel E. Bauer
bioRxiv 2024.07.24.605019; doi: https://doi.org/10.1101/2024.07.24.605019

## Features<a name="features"></a>

The package implements a two-stage pipeline comprising read preprocessing and off-target identification. The preprocessing module ingests pooled multi-sample FASTQ reads and demultiplexes them into sample-specific FASTQs for CRISPResso2 analysis.

<img src="./figure/absolveseq.svg" alt="absolveseq_flowchart" width="350">

## Dependencies<a name="dependencies"></a>
* Python 3
* [`CRISPResso2`](<https://github.com/pinellolab/CRISPResso2>) genome editing outcomes analysis tool

## Features<a name="features"></a>

The package implements a pipeline consisting of a read preprocessing module followed by an off-target identification module. The preprocessing module takes raw reads (FASTQ) from a pooled multi-sample sequencing run as input. Reads are demultiplexed into sample-specific FASTQs and plasmid barcode are moved to Header line.

The individual pipeline steps are:

1. **Target Barcode Demultiplexing:** Deconvolution of pooled sequencing data into target-specific read files.
2. **Plasmid Barcode Demultiplexing**: Stratification of target-specific reads by unique plasmid barcodes.
3. **Plasmid Barcode Filtering**: Identification and selection of high-quality plasmid barcodes from the pool.
4. **CRISPResso2 Analysis**: Quantification of editing outcomes using the CRISPResso2 pipeline.
5. **Deduplication & Recombination Filtering**: Removal of PCR duplicates and recombination artifacts.
6. **Indel Estimation & Power Analysis**: Statistical assessment of insertion/deletion frequencies and experimental power.
7. **Visualization & Reporting**: Generation of summary plots and comprehensive data reports.

## Running Analysis Steps Individually<a name="individual_steps"></a>

```bash
python absolveseq/1_demultiplex_by_targetBarcode.py \
  --output_dir ./test/demultiplexed_tBC_fastq \
  --target_info ./test/data/target_info/LVOT_oligo_pool.xlsx \
  --sample_info ./test/data/target_info/NovaSeq3_sample_info_example.csv \
  --n_processes 8
```

```bash
python absolveseq/2_demultiplex_by_plasmidBarcode.py \
  --fastq_dir ./test/demultiplexed_tBC_fastq \
  --fa_out_folder ./demultiplexed_pBC_fastq \
  --crispress_input_folder ./test/CRISPResso_input_files \
  --amplicon_fn ./test/data/target_info/OT_guide_amplicon_seq.csv \
  --n_processes 8
```


