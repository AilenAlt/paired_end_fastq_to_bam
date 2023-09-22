# Snakemake paired-end NGS Data Processing Workflow

This repository contains a Snakemake workflow for processing Next-Generation Sequencing (NGS) data. The workflow automates several data processing steps, including trimming, quality control, alignment, duplicate marking, and filtering.

## Prerequisites

Before running the workflow, ensure that you have the following software/tools installed on your system:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [BBMap/BBduk](https://sourceforge.net/projects/bbmap/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Samtools](http://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/)

Additionally, you will need reference genome/index files for your specific analysis.

## Usage

1. Clone this repository:

   ```shell
   git clone https://github.com/AilenAlt/single_end_fastq_to_bam
   cd your-repo
2. Create a config.yaml file to specify the paths to adapters and other configuration options.
   ```shell
   # config.yaml
    Adapters: path/to/adapters.fa
3. Edit the Snakefile to customize the workflow to your specific needs, such as reference genome/index paths and other parameters.
   
   In the `Snakefile`, you will find a section where a list of SRR identifiers is defined. These SRR identifiers represent specific NGS data samples that you       intend to process with this workflow.
   To use this workflow for your own NGS data analysis, you should replace the existing SRR list with your own list of SRR identifiers. The SRR identifiers you     provide should correspond to the NGS data samples you want to process.
   
4. Run the workflow using Snakemake:
    ```shell
    snakemake --cores <number_of_cores> --use-conda
    ```
    Replace <number_of_cores> with the desired number of CPU cores.

## Workflow Steps

- **Trimming**: Trims adapter sequences from raw FASTQ files using BBduk.
- **Quality Control (QC)**: Performs QC analysis using FastQC on the trimmed data.
- **Alignment**: Aligns the trimmed data to a reference genome using Bowtie2.
- **Sorting**: Sorts the resulting SAM files using Samtools.
- **Duplicate Marking**: Marks duplicate reads using Picard tools.
- **Filtering**: Filters the aligned data based on mapping quality, flag, and other criteria using Samtools.

<p align="center">
  <img src="https://github.com/AilenAlt/single_end_fastq_to_bam/blob/main/dag.png?raw=true" alt="Snakemake rules">
</p>

  ## Output
  
The workflow generates processed and filtered BAM files in the align directory. QC reports can be found in the qc directory.
