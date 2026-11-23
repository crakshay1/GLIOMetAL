#!/usr/bin/env Rscript

normalization <- function(DATA_DIR, OUTPUT_DIR) {

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
        "filtered_gene_counts.tsv"
    )

    TPM_FILE <- file.path(
        DATA_DIR,
        "filtered_gene_tpm.tsv"
    )

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

    # 2. INITIAL VALIDATION

    stopifnot(
        identical(colnames(counts), colnames(tpm))
    )

    stopifnot(
        identical(rownames(counts), rownames(tpm))
    )

    stopifnot(
        !anyNA(counts),
        !anyNA(tpm),
        all(is.finite(as.matrix(counts))),
        all(is.finite(as.matrix(tpm))),
        all(counts >= 0),
        all(tpm >= 0)
    )

    # 3. CREATE DGEList OBJECT

    y <- DGEList(counts = counts)

    # 4. TMM NORMALIZATION

    y <- calcNormFactors(
        y,
        method = "TMM"
    )
    cat("\nTMM normalization factors:\n")
    print(y$samples)
    
    cat("\nGeometric mean of TMM factors:\n")
    print(exp(mean(log(y$samples$norm.factors))))
    
    cat("\nEffective library sizes:\n")
    print(y$samples$lib.size * y$samples$norm.factors)

    normalization_factors <- data.frame(
        Sample = colnames(y),
        Library_Size = y$samples$lib.size,
        Norm_Factor = y$samples$norm.factors,
        Effective_Library_Size =
            y$samples$lib.size * y$samples$norm.factors
    )

    # 5. LOG2-CPM TRANSFORMATION

    logCPM <- cpm(
        y,
        log = TRUE,
        prior.count = 2
    )

    # 6. TPM TRANSFORMATION

    logTPM <- log2(
        tpm + 1
    )

    # 7. EXPORT RESULTS

    write.table(
        logCPM,
        file.path(OUTPUT_DIR, "normalized_log2CPM.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    write.table(
        logTPM,
        file.path(OUTPUT_DIR, "log2TPM.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    write.table(
        normalization_factors,
        file.path(OUTPUT_DIR, "normalization_factors.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("\nNormalization completed successfully!\n")
    cat("Normalization method: TMM\n")
    cat("Transformation: log2-CPM\n")
    cat("Additional transformation: log2(TPM + 1)\n")
}



args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop(
        "Usage: Rscript normalization.R <DATA_DIR> <OUTPUT_DIR>"
    )
}

normalization(
    DATA_DIR = args[1],
    OUTPUT_DIR = args[2]
)
