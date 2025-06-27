# ParaMAG: Pipeline for Metagenome-Assembled Genomes for Intracellular Parasites

This repository contains a Snakemake pipeline and a corresponding conda environment file for performing metagenomic assembly, binning, and phylogenomic analysis, with a focus on intracellular parasites and other organisms with reduced genomes.

## Requirements
- Snakemake
- Python 3.8+
- Conda (optional, but recommended)

## Installation
Clone this repository:
```
git clone https://github.com/mikeraiko/ParaMAG.git
cd ParaMAG
```

## Usage

### Step 1: Run with your own data
Prepare your input files and create a minimal config file `config.yaml` like this:
```yaml
reads1: data/sample_R1.fastq.gz
reads2: data/sample_R2.fastq.gz
adapter_file: adapters/TruSeq3-PE.fa
```
Then run the wrapper script:
```
python paramag.py --run --config config.yaml
```



### Step 2: Run on example data
To test the pipeline:
```
python paramag.py --example
```
This will automatically run the full workflow on a small test dataset.

### Help
For help and additional options:
```
python paramag.py --help
```

## Output
- `maxbin_output/`: binned genomes
- `busco_output/`: completeness metrics
- `super_tree_output.tre`: consensus tree
- `report/summary_report.md`: summary report
