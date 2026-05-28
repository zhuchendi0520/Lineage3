# Lineage 3 Analysis Scripts

This repository contains scripts and supporting files used for the analysis of *Mycobacterium tuberculosis* Lineage 3 (L3), including variant calling, sublineage assignment, population structure analysis, phylogenetic reconstruction, evolutionary dating, demographic history, terminal branch length analysis, transmission cluster inference, drug-resistance mutation screening, pNS calculation, and THD analysis.

## Repository Structure

| Directory | Description |
| --- | --- |
| `Variant Calling and Genome Assembly/` | Scripts and reference files for read trimming, mapping, SNP calling, filtering, and SNP-format conversion. |
| `L3 sublineage assign/` | Perl script for assigning L3.1.1 sublineages using diagnostic SNP positions. |
| `Population Structure and Phylogenetic Reconstruction/` | R script for fastBAPS population structure analysis and visualization with phylogenetic trees. |
| `Evolutionary Dating and Demographic History/` | BEAST XML files and R scripts for evolutionary dating and Bayesian skyline plot analysis. |
| `Phylogeographic Reconstruction/` | Phylogeographic reconstruction output files and compressed archives for L3.1.1 sublineages. |
| `Phylogeographic Inference of Migration-Associated Transmission/` | Python script for identifying transmission clusters and updating migration-associated classification from tree and metadata files. |
| `Terminal Branch Length (TBL) Analysis/` | Scripts for extracting terminal branch lengths from tree files and downstream TBL analysis. |
| `Genotypic Resistance Prediction and Compensatory Mutation Analysis/` | Scripts and mutation database files for identifying drug-resistance mutations and compensatory mutations from SNP files. |
| `pNS/` | Python scripts for calculating pNS values from annotation files. |
| `THD/` | R script for calculating transmission fitness / THD-related indices and generating plots. |

## Requirements

The scripts use a mixture of Python, R, Perl, and external bioinformatics tools.

### Python

Recommended Python version: Python 3.

Common Python packages:

```bash
pip install biopython
```

Some scripts only use the Python standard library.

### R

Common R packages used by the analysis scripts include:

```r
install.packages(c(
  "ape",
  "ggplot2",
  "data.table",
  "stringr",
  "readxl",
  "lmerTest",
  "beeswarm",
  "devtools"
))

devtools::install_github("gtonkinhill/fastbaps")
```

The THD analysis additionally requires the `thd` package.

### External Tools

Variant calling scripts are written to generate or run workflows using tools such as:

- sickle
- bowtie2
- samtools
- VarScan
- Perl

Please update tool paths and reference genome paths in the scripts before running them in a new computing environment.

## Usage Examples

### 1. L3.1.1 Sublineage Assignment

Assign an L3.1.1 sublineage based on diagnostic SNP positions:

```bash
perl "L3 sublineage assign/lineage3.1.pl" sample.snp
```

Diagnostic SNPs used in the script:

| Position | Assigned sublineage |
| --- | --- |
| 2118846 | L3.1.1.1 |
| 1914217 | L3.1.1.2 |
| 2812520 | L3.1.1.3 |
| 3350236 | L3.1.1.4 |
| 1134143 | L3.1.1.5 |
| 3592529 | L3.1.1.6 |

### 2. Variant Calling and Genome Assembly

Generate mapping and SNP-calling shell commands for paired-end sequencing data:

```bash
python "Variant Calling and Genome Assembly/1_mapping_pair_fixed_nostrandbias.py" strain_list.txt
```

Generate mapping and SNP-calling shell commands for single-end sequencing data:

```bash
python "Variant Calling and Genome Assembly/1_mapping_single_fixed_nostrandbias.py" strain_list.txt
```

The input `strain_list.txt` should contain one strain/sample ID per line. The scripts expect FASTQ files named according to the sample ID.

### 3. Population Structure Analysis

Run fastBAPS-based clustering and tree visualization:

```r
source("Population Structure and Phylogenetic Reconstruction/fastBAPS.R")
```

Before running, update the working directory, input FASTA file, and tree file names in the R script.

### 4. Terminal Branch Length Extraction

Extract terminal branch lengths from a Newick tree file:

```bash
python "Terminal Branch Length (TBL) Analysis/extract_tree_TBL.py" input.tree
```

The script writes:

```text
input_TBL.csv
```

### 5. Migration-Associated Transmission Cluster Inference

Infer clusters from a Newick tree and metadata table:

```bash
python "Phylogeographic Inference of Migration-Associated Transmission/infer_transmission_clusters.py" tree.nwk 0.0005 metadata.csv
```

The metadata file should contain at least the following columns:

```text
Strain,Country,Migration
```

The output file is named:

```text
tree_cluster_migration.csv
```

### 6. Drug-Resistance Mutation Extraction

Extract drug-resistance mutations from a mutation database and SNP file:

```bash
python "Genotypic Resistance Prediction and Compensatory Mutation Analysis/snpfile2extract_dr_mutation.py" \
  "Genotypic Resistance Prediction and Compensatory Mutation Analysis/Total_drug_resistance_mutation.txt" \
  sample.snp
```

### 7. pNS Analysis

Calculate pNS values from annotation CSV files:

```bash
python pNS/pNS.py sample.ann.csv
```

If a gene has no synonymous mutation and the synonymous count should be set to 1, use:

```bash
python pNS/pNS_setsynto1.py sample.ann.csv
```

### 8. THD Analysis

Run the THD analysis script in R:

```r
source("THD/THD_L3.R")
```

Before running, update the working directory and input file names in the script.

## Notes

- Several scripts contain absolute paths from the original analysis environment. These should be replaced with paths appropriate for the current system.
- Input file names, metadata columns, and FASTQ naming conventions should be checked before running each module.
- Large intermediate files and raw sequencing data are not included in this repository.
- BEAST XML files and phylogeographic output files are provided as analysis configuration/results files for reproducibility.

## Citation

If these scripts are used in a publication or shared analysis workflow, please cite the associated study or repository.
