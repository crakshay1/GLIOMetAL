#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def quant_qc(DATA_DIR):
    """
        This module gives the library size, the detected genes, the zero count proportion & the TPM distribution.
    """
    # CONFIGURATION
    SALMON_DIR = DATA_DIR / "Data_for_M2_Biomede"

    OUTPUT_DIR = Path("QUANTIFICATION")
    OUTPUT_DIR.mkdir(exist_ok=True)

    GENE_COUNTS = DATA_DIR / "merged_samples_SYMBOL_raw_count.tsv"
    GENE_TPM = DATA_DIR / "merged_samples_SYMBOL_tpm.tsv"

    TRANSCRIPT_COUNTS = DATA_DIR / "merged_samples_txID_raw_count.tsv"
    TRANSCRIPT_TPM = DATA_DIR / "merged_samples_txID_tpm.tsv"


    # LOAD MATRICES
    gene_counts = pd.read_csv(GENE_COUNTS, sep="\t", index_col=0)
    gene_tpm = pd.read_csv( GENE_TPM, sep="\t", index_col=0)
    transcript_counts = pd.read_csv(TRANSCRIPT_COUNTS, sep="\t", index_col=0)
    transcript_tpm = pd.read_csv(TRANSCRIPT_TPM, sep="\t", index_col=0)


    # SAMPLE-LEVEL QUANTIFICATION QC
    qc = pd.DataFrame(index=gene_counts.columns)

    # Library size
    qc["total_gene_counts"] = gene_counts.sum(axis=0)
    qc["total_transcript_counts"] = transcript_counts.sum(axis=0)

    # Number of detected features
    qc["detected_genes"] = (gene_counts > 0).sum(axis=0)
    qc["detected_transcripts"] = (
        transcript_counts > 0
    ).sum(axis=0)


    # Zero-count proportion
    qc["gene_zero_proportion"] = (
        (gene_counts == 0).sum(axis=0) / len(gene_counts)
    )
    qc["transcript_zero_proportion"] = (
        (transcript_counts == 0).sum(axis=0)
        / len(transcript_counts)
    )


    # Total TPM
    qc["gene_total_tpm"] = gene_tpm.sum(axis=0)
    qc["transcript_total_tpm"] = transcript_tpm.sum(axis=0)


    # TPM DISTRIBUTION STATISTICS
    tpm_stats = pd.DataFrame(index=gene_tpm.columns)
    tpm_stats["mean_TPM"] = gene_tpm.mean(axis=0)
    tpm_stats["median_TPM"] = gene_tpm.median(axis=0)
    tpm_stats["Q1_TPM"] = gene_tpm.quantile(0.25, axis=0)
    tpm_stats["Q3_TPM"] = gene_tpm.quantile(0.75, axis=0)
    tpm_stats["max_TPM"] = gene_tpm.max(axis=0)
    tpm_stats["mean_log2_TPM"] = np.log2(
        gene_tpm + 1
    ).mean(axis=0)


    # SAVE RESULTS
    qc.to_csv(OUTPUT_DIR / "sample_qc_metrics.tsv", sep="\t")
    tpm_stats.to_csv(
        OUTPUT_DIR / "tpm_distribution_statistics.tsv",
        sep="\t"
    )

    print("\n===== SAMPLE QC METRICS =====")
    print(qc.round(4))
    print("\n===== TPM DISTRIBUTION =====")
    print(tpm_stats.round(4))


    # VISUALIZATION 1: LIBRARY SIZE
    qc["total_gene_counts"].sort_values().plot(
        kind="bar",
        figsize=(12, 5)
    )

    plt.ylabel("Total estimated counts")
    plt.xlabel("Sample")
    plt.title("Library size per sample")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "library_size.png", dpi=300)
    plt.close()


    # VISUALIZATION 2: DETECTED GENES
    qc["detected_genes"].sort_values().plot(
        kind="bar",
        figsize=(12, 5)
    )

    plt.ylabel("Number of detected genes")
    plt.xlabel("Sample")
    plt.title("Detected genes per sample")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "detected_genes.png", dpi=300)
    plt.close()


    # VISUALIZATION 3: ZERO-COUNT PROPORTION
    qc["gene_zero_proportion"].sort_values().plot(
        kind="bar",
        figsize=(12, 5)
    )

    plt.ylabel("Proportion of zero counts")
    plt.xlabel("Sample")
    plt.title("Gene-level zero-count proportion")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "zero_count_proportion.png", dpi=300)
    plt.close()


    # VISUALIZATION 4: TPM DISTRIBUTION
    plt.figure(figsize=(12, 6))
    for sample in gene_tpm.columns:
        values = np.log2(gene_tpm[sample] + 1)
        plt.hist(
            values,
            bins=60,
            density=True,
            histtype="step",
            linewidth=1,
            label=sample
        )

    plt.xlabel("log2(TPM + 1)")
    plt.ylabel("Density")
    plt.title("Gene expression distributions across samples")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "tpm_distributions.png", dpi=300)
    plt.close()


    print(f"\nQC results saved in: {OUTPUT_DIR}")
