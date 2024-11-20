# ParaMAG: Pipeline for Metagenome-Assembled Genomes for Intracellular Parasites

This repository contains a Snakemake pipeline and a corresponding conda environment file for performing metagenomic assembly, binning, and phylogenomic analysis, with a focus on intracellular parasites and other organisms with reduced genomes.



### Pipeline Overview
Steps:
1. **Quality Control:**
Uses FastQC to evaluate read quality and Trimmomatic for adapter removal and trimming.
2. **Assembly:**
Assembles trimmed reads into scaffolds using SPAdes in single-cell mode.
3. **Binning:**
Groups assembled scaffolds into bins using MaxBin2.
4. **Completeness Check:**
Evaluates the completeness of each bin with BUSCO.
5. **Phylogenomic Analysis:**
Aligns sequences with MAFFT, trims alignments with trimAL, and performs phylogenetic tree inference with IQ-TREE.
Optionally, reconstructs supertrees with Astral using gene trees.
### Requirements
The following tools are required for this pipeline:
- FastQC: Read quality control
- Trimmomatic: Read trimming
- SPAdes: Genome assembly
- MaxBin2: Genome binning
- BUSCO: Bin completeness assessment
- MAFFT: Multiple sequence alignment
- trimAL: Alignment trimming
- IQ-TREE: Phylogenetic inference
- Astral: Supertree reconstruction (optional)
- Python 3.6 or higher
- Snakemake (version ≥ 5.0)


### Setup Instructions
1. Clone this repository to your local machine:
```
git clone <repository_url>
cd <repository_folder>
```
2. Install Miniconda or Conda (if not already installed). You can download it from [here](https://docs.anaconda.com/miniconda/).

3. Use the provided environment.yaml file to create the conda environment:

`conda env create -f environment.yaml`

4. Activate the environment to use the installed dependencies:

`conda activate paramag`

### Running the Pipeline
1. Prepare Input Files
 
Place your paired-end raw sequencing reads in the input_reads/ directory:

File names: sample_R1.fastq and sample_R2.fastq.

Ensure the fungi_odb10 database for BUSCO is downloaded and accessible.

2. Run Snakemake

To execute the pipeline:

`snakemake -p --cores <number_of_cores>`

The -p flag prints each command before execution.

Replace <number_of_cores> with the number of parallel threads you want to use.

3. View Output
  
Key outputs include:

- Bins: maxbin_output/maxbin_output.*.fasta.
- Completeness Reports: busco_output/busco_{n}/full_table.tsv.
- Taxonomy Results: taxonomy_results/bin_{n}_taxonomy.txt.
- Phylogenomic Trees: Trees generated from alignment and trimming steps.

### Reproducibility
To ensure reproducibility, use the same environment.yaml file to recreate the conda environment.

### License
This pipeline is open-source and distributed under the MIT License.

### Support
For questions or bug reports, please open an issue in this repository or contact me at [mike.rayko@gmail.com].
