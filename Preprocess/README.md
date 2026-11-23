# Preprocessing


**Quality control, filtering, normalization and descriptive statistics of bulk RNA-seq data (Salmon quantification)**  
*Digital Genomics — Université Paris-Saclay*

-----------------------------------------------------

This pipeline prepares an RNA-seq dataset for downstream differential analysis. It chains Python modules (metadata, integrity and quantification checks) and R scripts (filtering, normalization, descriptive statistics) from a single entry point, `preprocess.py`. The experimental groups are defined by the `TP53_status` variable (two groups, 15 samples each, 30 samples in total).

## Input data

The input folder is selected through a file dialog (easygui) when the pipeline starts. It must contain:

```
Metadata_for_M2_Biomede.txt                 # Sample metadata (tab-separated)
merged_samples_SYMBOL_raw_count.tsv         # Gene-level raw counts
merged_samples_SYMBOL_tpm.tsv               # Gene-level TPM
merged_samples_txID_raw_count.tsv           # Transcript-level raw counts
merged_samples_txID_tpm.tsv                 # Transcript-level TPM
Data_for_M2_Biomede/<sample>/quant.sf       # Salmon outputs
```

## Pipeline

| Step | Script | Description | Output folder |
|------|--------|-------------|---------------|
| 1 | `Metadata.py` | Metadata QC: missing values, constant and duplicated variables, categorical and numerical summaries, consistency with the expression matrix | `metadata_qc/` |
| 2 | `Integrity.py` | Integrity of the four matrices and of the Salmon `quant.sf` files: missing, infinite, negative and non-numeric values, duplicated identifiers, sample consistency | `INTEGRITY/` |
| 3 | `Quantification.py` | Library size, detected genes, zero-count proportion and TPM distributions per sample | `QUANTIFICATION/` |
| 4 | `Filter.r` | Low-expression filtering with edgeR `filterByExpr`, matching metadata with samples and checking the TP53 groups | `FILTERING/` |
| 5 | `Normalization.r` | TMM normalization, log2-CPM and log2(TPM + 1) transformations | `NORMALIZATION/` |
| 6 | `Desc.r` | Sample and gene statistics, top variable genes, sample correlation, PCA and expression distributions | `DESCRIPTIVE/` |

The R scripts are called from `preprocess.py` with `os.system`. Before step 6, the TPM matrix and the metadata file are copied into `NORMALIZATION/`, because `Desc.r` reads all its inputs from a single folder.

## Tech stack

* Python (pandas, numpy, matplotlib, easygui)
* R (edgeR, ggplot2, pheatmap)
* Salmon quantification files as input

## Installation

Python dependencies:

```bash
pip install -r requirements.txt
```

R packages (the system libraries are needed to compile edgeR on Debian/Ubuntu):

```bash
sudo apt install -y build-essential gfortran libblas-dev liblapack-dev
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages(c("ggplot2","pheatmap","BiocManager")); BiocManager::install("edgeR", ask=FALSE, update=FALSE)'
```

## Usage

```bash
python preprocess.py
```

A dialog opens to select the input folder. All output folders are created in the current working directory. Each R script can also be run on its own:

```bash
Rscript Filter.r <DATA_DIR> <OUTPUT_DIR>
Rscript Normalization.r <FILTERED_DIR> <OUTPUT_DIR>
Rscript Desc.r <NORMALIZED_DIR> <OUTPUT_DIR>
```

## Filtering and normalization parameters

* Filtering: `filterByExpr` with `min.count = 10`, `min.total.count = 15` and `min.n = 3`, using the TP53 groups
* Normalization: TMM factors, log2-CPM with `prior.count = 2`
* Descriptive statistics: PCA on the 500 most variable genes

## Structure

```
.
├── preprocess.py                   # Entry point
├── Metadata.py                     # Metadata QC
├── Integrity.py                    # Matrix and Salmon integrity checks
├── Quantification.py               # Quantification QC
├── Filter.r                        # Low-expression filtering
├── Normalization.r                 # TMM normalization and transformations
├── Desc.r                          # Descriptive statistics and figures
├── metadata_qc/                    # Metadata QC tables and figures
├── INTEGRITY/
│	 └── integrity_report.csv        # PASS / WARNING / FAIL report
├── QUANTIFICATION/
│	├── sample_qc_metrics.tsv
│	├── tpm_distribution_statistics.tsv
│	├── library_size.png
│	├── detected_genes.png
│	├── zero_count_proportion.png
│	└── tpm_distributions.png
├── FILTERING/
│	├── filtered_gene_counts.tsv
│	├── filtered_gene_tpm.tsv
│	├── gene_filtering_status.tsv
│	└── constant_genes.tsv
├── NORMALIZATION/
│	├── normalized_log2CPM.tsv
│	├── log2TPM.tsv
│	└── normalization_factors.tsv
├── DESCRIPTIVE/
│	├── sample_descriptive_stats.tsv
│	├── sample_statistics_by_group.tsv
│	├── gene_descriptive_stats.tsv
│	├── top_variable_genes.tsv
│	├── sample_correlation.tsv
│	├── PCA_coordinates.tsv
│	└── PCA_explained_variance.tsv
├── requirements.txt            # Python dependencies
└── figures/                    # QC, PCA, correlation, distributions, group comparisons
```
