#!/usr/bin/env Rscript

filtering <- function(DATA_DIR, OUTPUT_DIR) {

    suppressMessages(
        suppressWarnings(library(edgeR))
    )

    dir.create(
        OUTPUT_DIR,
        showWarnings = FALSE,
        recursive = TRUE
    )

    COUNTS_FILE <- file.path(
        DATA_DIR,
        "merged_samples_SYMBOL_raw_count.tsv"
    )

    TPM_FILE <- file.path(
        DATA_DIR,
        "merged_samples_SYMBOL_tpm.tsv"
    )

    METADATA_FILE <- file.path(
        DATA_DIR,
        "Metadata_for_M2_Biomede.txt"
    )

    SAMPLE_ID_COLUMN <- "sample"
    GROUP_COLUMN <- "TP53_status"

    EXPECTED_SAMPLES <- 30

    MIN_SAMPLES <- 3
    MIN_COUNT <- 10
    MIN_TOTAL_COUNT <- 15

    # 1. LOAD DATA

    counts <- read.delim(
        COUNTS_FILE,
        row.names = 1,
        check.names = FALSE
    )

    tpm <- read.delim(
        TPM_FILE,
        row.names = 1,
        check.names = FALSE
    )

    metadata <- read.delim(
        METADATA_FILE,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    # 2. INITIAL VALIDATION

    stopifnot(
        identical(colnames(counts), colnames(tpm))
    )

    stopifnot(
        !anyDuplicated(rownames(counts)),
        !anyDuplicated(rownames(tpm))
    )

    stopifnot(
        setequal(rownames(counts), rownames(tpm))
    )

    if (!SAMPLE_ID_COLUMN %in% colnames(metadata)) {
        stop("Sample ID column not found in metadata.")
    }

    if (!GROUP_COLUMN %in% colnames(metadata)) {
        stop("Group column not found in metadata.")
    }

    if (anyDuplicated(metadata[[SAMPLE_ID_COLUMN]])) {
        stop("Duplicated sample identifiers in metadata.")
    }

    tpm <- tpm[rownames(counts), , drop = FALSE]

    if (anyNA(counts) || any(!is.finite(as.matrix(counts)))) {
        stop("Invalid values detected in count matrix.")
    }

    if (any(counts < 0)) {
        stop("Negative counts detected.")
    }

    if (anyNA(tpm) || any(!is.finite(as.matrix(tpm)))) {
        stop("Invalid values detected in TPM matrix.")
    }

    if (any(tpm < 0)) {
        stop("Negative TPM values detected.")
    }

    # 3. MATCH METADATA WITH EXPRESSION MATRIX

    expression_samples <- colnames(counts)

    metadata_index <- match(
        expression_samples,
        metadata[[SAMPLE_ID_COLUMN]]
    )

    if (anyNA(metadata_index)) {
        stop("Some expression samples are missing from metadata.")
    }

    metadata_ordered <- metadata[
        metadata_index,
        ,
        drop = FALSE
    ]

    if (ncol(counts) != EXPECTED_SAMPLES) {
        warning(
            "Unexpected number of expression samples: ",
            ncol(counts)
        )
    }

    GROUP <- factor(
        metadata_ordered[[GROUP_COLUMN]]
    )

    if (anyNA(GROUP)) {
        stop("Missing group assignments detected.")
    }

    if (nlevels(GROUP) != 2) {
        stop("Expected exactly two TP53 groups.")
    }

    cat("\n===== EXPERIMENTAL GROUPS =====\n")
    print(table(GROUP))

    if (!identical(
        sort(as.integer(table(GROUP))),
        c(15L, 15L)
    )) {
        warning(
            "Group sizes differ from the expected 15/15 distribution."
        )
    }

    # 4. INITIAL STATISTICS

    initial_genes <- nrow(counts)

    cat("Initial number of genes:", initial_genes, "\n")

    all_zero <- rowSums(counts) == 0

    cat("All-zero genes:", sum(all_zero), "\n")

    # 5. CREATE DGEList OBJECT

    y <- DGEList(
        counts = counts,
        group = GROUP
    )

    # 6. LOW-EXPRESSION FILTERING

    keep <- filterByExpr(
        y,
        group = GROUP,
        min.count = MIN_COUNT,
        min.total.count = MIN_TOTAL_COUNT,
        min.n = MIN_SAMPLES
    )

    filtered_counts <- counts[
        keep,
        ,
        drop = FALSE
    ]

    filtered_tpm <- tpm[
        keep,
        ,
        drop = FALSE
    ]

    # 7. FILTERING SUMMARY

    final_genes <- nrow(filtered_counts)
    removed_genes <- initial_genes - final_genes

    retained_percentage <- (
        final_genes / initial_genes
    ) * 100

    cat("\n===== FILTERING SUMMARY =====\n")

    cat("Initial genes:", initial_genes, "\n")
    cat("Retained genes:", final_genes, "\n")
    cat("Removed genes:", removed_genes, "\n")

    cat(
        "Retained percentage:",
        round(retained_percentage, 2),
        "%\n"
    )

    # 8. IDENTIFY CONSTANT GENES

    gene_variance <- apply(
        filtered_counts,
        1,
        var
    )

    constant_genes <- gene_variance == 0

    cat(
        "Constant genes after filtering:",
        sum(constant_genes),
        "\n"
    )

    # 9. EXPORT RESULTS

    write.table(
        filtered_counts,
        file.path(OUTPUT_DIR, "filtered_gene_counts.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    write.table(
        filtered_tpm,
        file.path(OUTPUT_DIR, "filtered_gene_tpm.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    write.table(
        data.frame(
            Gene = rownames(counts),
            Keep = keep
        ),
        file.path(OUTPUT_DIR, "gene_filtering_status.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    write.table(
        data.frame(
            Gene = rownames(filtered_counts)[constant_genes],
            Variance = gene_variance[constant_genes]
        ),
        file.path(OUTPUT_DIR, "constant_genes.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("\nFiltering completed successfully!\n")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript filtering.R <DATA_DIR> <OUTPUT_DIR>")
}

filtering(
    DATA_DIR = args[1],
    OUTPUT_DIR = args[2]
)
