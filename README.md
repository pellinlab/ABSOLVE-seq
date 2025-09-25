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

1. **Sample demultiplexing and plasmid barcode extraction**: A pooled multi-sample sequencing run is demultiplexed into sample-specific read files based on sample-specific
2. **Aggragate reads with plasmid barcode**:
3. **CRISPResso2 analyais**:
4. **Dedud and recombination filtering**: 
5. **Indel estimation and power analysis**:
6. **Reporting**:
7. **Visualization**:


## Running Analysis Steps Individually<a name="individual_steps"></a>