#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
})

descriptive_stats <- function(DATA_DIR, OUTPUT_DIR) {

    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

    FIG <- file.path(OUTPUT_DIR, "figures")

	for (folder in c("expression", "samples", "clinical",
	                 "molecular", "groups")) {
	    dir.create(file.path(FIG, folder), recursive = TRUE, showWarnings = FALSE)
	}

    LOGCPM_FILE <- file.path(DATA_DIR, "normalized_log2CPM.tsv")
    TPM_FILE <- file.path(DATA_DIR, "merged_samples_SYMBOL_tpm.tsv")
    METADATA_FILE <- file.path(DATA_DIR, "Metadata_for_M2_Biomede.txt")

    SAMPLE_ID <- "sample"
    GROUP_COLUMN <- "DMG_subgroup"
    TOP_GENES <- 500

    # 1. LOAD DATA

    logCPM <- read.delim(LOGCPM_FILE, row.names = 1, check.names = FALSE)
    tpm <- read.delim(TPM_FILE, row.names = 1, check.names = FALSE)
    metadata <- read.delim(METADATA_FILE, check.names = FALSE)

    logCPM <- as.matrix(logCPM)
    tpm <- as.matrix(tpm)

    storage.mode(logCPM) <- "numeric"
    storage.mode(tpm) <- "numeric"

    if (anyDuplicated(rownames(logCPM)) ||
        anyDuplicated(rownames(tpm))) {
        stop("Duplicated gene identifiers.")
    }

    if (!all(colnames(logCPM) %in% colnames(tpm))) {
        stop("Sample mismatch between logCPM and TPM.")
    }

    if (!all(rownames(logCPM) %in% rownames(tpm))) {
        stop("Some logCPM genes are missing from TPM.")
    }

    if (!all(colnames(logCPM) %in% metadata[[SAMPLE_ID]])) {
        stop("Expression samples missing from metadata.")
    }

    metadata <- metadata[
        match(colnames(logCPM), metadata[[SAMPLE_ID]]), ,
        drop = FALSE
    ]

    if (anyNA(metadata[[GROUP_COLUMN]])) {
        stop("Missing DMG_subgroup information.")
    }

    GROUP <- factor(metadata[[GROUP_COLUMN]])

    if (nlevels(GROUP) != 2) {
        warning("The number of DMG subgroups is not equal to two.")
    }

    tpm <- tpm[rownames(logCPM), colnames(logCPM), drop = FALSE]

    # 2. SAMPLE STATISTICS

    detected <- colSums(logCPM > 1)

    sample_stats <- data.frame(
        Sample = colnames(logCPM),
        Group = GROUP,
        Mean = colMeans(logCPM),
        Median = apply(logCPM, 2, median),
        SD = apply(logCPM, 2, sd),
        Q1 = apply(logCPM, 2, quantile, 0.25),
        Q3 = apply(logCPM, 2, quantile, 0.75),
        Detected_Genes = detected,
        Detected_Genes_Percent = detected / nrow(logCPM) * 100
    )

    write.table(
        sample_stats,
        file.path(OUTPUT_DIR, "sample_descriptive_stats.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    group_stats <- aggregate(
        sample_stats[, c("Mean", "Median", "SD",
                         "Detected_Genes", "Detected_Genes_Percent")],
        list(Group = GROUP),
        function(x) c(mean = mean(x), median = median(x), sd = sd(x))
    )

    write.table(
        group_stats,
        file.path(OUTPUT_DIR, "sample_statistics_by_group.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # 3. GENE STATISTICS

    gene_stats <- data.frame(
        Gene = rownames(logCPM),
        Mean_logCPM = rowMeans(logCPM),
        Median_logCPM = apply(logCPM, 1, median),
        SD_logCPM = apply(logCPM, 1, sd),
        Variance_logCPM = apply(logCPM, 1, var),
        MAD_logCPM = apply(logCPM, 1, mad),
        Mean_TPM = rowMeans(tpm),
        SD_TPM = apply(tpm, 1, sd)
    )

    gene_stats$CV_TPM <- ifelse(
        gene_stats$Mean_TPM > 0,
        gene_stats$SD_TPM / gene_stats$Mean_TPM,
        NA_real_
    )

    gene_stats <- gene_stats[
        order(gene_stats$Variance_logCPM, decreasing = TRUE), ]

    write.table(
        gene_stats,
        file.path(OUTPUT_DIR, "gene_descriptive_stats.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
        head(gene_stats, TOP_GENES),
        file.path(OUTPUT_DIR, "top_variable_genes.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # 4. SAMPLE CORRELATION

    sample_cor <- cor(
        logCPM,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    write.table(
        sample_cor,
        file.path(OUTPUT_DIR, "sample_correlation.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    annotation <- data.frame(Group = GROUP)
    rownames(annotation) <- colnames(logCPM)

    heatmap <- pheatmap(
        sample_cor,
        annotation_col = annotation,
        annotation_row = annotation,
        main = "Sample-to-sample Pearson correlation",
        silent = TRUE
    )

    png(
        file.path(FIG, "samples", "sample_correlation.png"),
        width = 1800,
        height = 1600,
        res = 200
    )

    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)

    dev.off()


    # 5. PCA

    gene_variances <- apply(logCPM, 1, var)

    selected_genes <- names(
        sort(gene_variances, decreasing = TRUE)
    )[seq_len(min(TOP_GENES, length(gene_variances)))]

    pca <- prcomp(
        t(logCPM[selected_genes, , drop = FALSE]),
        center = TRUE,
        scale. = FALSE
    )

    explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

    pca_coordinates <- data.frame(
        Sample = rownames(pca$x),
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        Group = GROUP,
        metadata[, setdiff(
            names(metadata),
            c(SAMPLE_ID, "files", "num", GROUP_COLUMN)
        ), drop = FALSE]
    )

    write.table(
        pca_coordinates,
        file.path(OUTPUT_DIR, "PCA_coordinates.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
        data.frame(
            Component = paste0("PC", seq_along(explained)),
            Explained_Variance_Percent = explained
        ),
        file.path(OUTPUT_DIR, "PCA_explained_variance.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    pca_plot <- ggplot(
        pca_coordinates,
        aes(PC1, PC2, color = Group)
    ) +
        geom_point(size = 3) +
        labs(
            title = "Principal Component Analysis",
            x = paste0("PC1 (", round(explained[1], 2), "%)"),
            y = paste0("PC2 (", round(explained[2], 2), "%)")
        ) +
        theme_bw()

    ggsave(
        file.path(FIG, "samples", "PCA.png"),
        pca_plot, width = 8, height = 6, dpi = 300
    )

    # 6. PCA BY METADATA

    pca_variables <- intersect(
        c("Sexe", "kit", "Localisation", "data_omics"),
        names(metadata)
    )

    if ("age" %in% names(metadata)) {
        pca_coordinates$Age_category <- cut(
            as.numeric(metadata$age),
            breaks = 3,
            labels = c("Low", "Intermediate", "High")
        )
        pca_variables <- c(pca_variables, "Age_category")
    }

    for (variable in pca_variables) {

        pca_coordinates$Color <- as.factor(
            pca_coordinates[[variable]]
        )

        plot <- ggplot(
            pca_coordinates,
            aes(PC1, PC2, color = Color)
        ) +
            geom_point(size = 3) +
            labs(
                title = paste("PCA colored by", variable),
                x = paste0("PC1 (", round(explained[1], 2), "%)"),
                y = paste0("PC2 (", round(explained[2], 2), "%)"),
                color = variable
            ) +
            theme_bw()

        ggsave(
            file.path(FIG, "samples", paste0("PCA_", variable, ".png")),
            plot, width = 8, height = 6, dpi = 300
        )
    }

    # 7. DENDROGRAM

    hc <- hclust(as.dist(1 - sample_cor), method = "average")

    png(file.path(FIG, "samples", "sample_dendrogram.png"),
        width = 1800, height = 1200, res = 200)

    plot(
        hc,
        main = "Sample hierarchical clustering",
        xlab = "",
        sub = ""
    )

    dev.off()

    # 8. EXPRESSION DISTRIBUTIONS

    expression_df <- data.frame(
        Expression = as.vector(logCPM),
        Sample = rep(colnames(logCPM), each = nrow(logCPM)),
        Group = rep(GROUP, each = nrow(logCPM))
    )

    distribution_plot <- ggplot(
        expression_df,
        aes(Sample, Expression, fill = Group)
    ) +
        geom_boxplot(outlier.size = 0.1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(
            title = "Expression distribution across samples",
            y = "log2-CPM",
            x = "Sample"
        )

    ggsave(
        file.path(FIG, "expression", "expression_distribution.png"),
        distribution_plot, width = 12, height = 6, dpi = 300
    )

    density_plot <- ggplot(
        expression_df,
        aes(Expression, group = Sample, color = Group)
    ) +
        geom_density(alpha = 0.15) +
        theme_bw() +
        labs(
            title = "Expression density across samples",
            x = "log2-CPM",
            y = "Density"
        )

    ggsave(
        file.path(FIG, "expression", "expression_density.png"),
        density_plot, width = 8, height = 6, dpi = 300
    )

    group_expression <- ggplot(
        expression_df,
        aes(Group, Expression, fill = Group)
    ) +
        geom_violin(trim = TRUE) +
        geom_boxplot(width = 0.12, outlier.size = 0.1) +
        theme_bw() +
        labs(
            title = "Expression distribution by DMG subgroup",
            x = "DMG subgroup",
            y = "log2-CPM"
        )

    ggsave(
        file.path(FIG, "groups", "expression_by_group.png"),
        group_expression, width = 7, height = 6, dpi = 300
    )

    # 9. MEAN-VARIANCE PLOT

    mean_variance <- ggplot(
        gene_stats,
        aes(Mean_logCPM, Variance_logCPM)
    ) +
        geom_point(alpha = 0.3, size = 0.7) +
        scale_y_log10() +
        theme_bw() +
        labs(
            title = "Mean-variance relationship",
            x = "Mean log2-CPM",
            y = "Variance of log2-CPM"
        )

    ggsave(
        file.path(FIG, "expression", "mean_variance.png"),
        mean_variance, width = 8, height = 6, dpi = 300
    )

    # 10. DETECTED GENES

    detection_plot <- ggplot(
        sample_stats,
        aes(Group, Detected_Genes, fill = Group)
    ) +
        geom_boxplot() +
        geom_jitter(width = 0.1, size = 2) +
        theme_bw() +
        labs(
            title = "Detected genes by DMG subgroup",
            x = "DMG subgroup",
            y = "Number of detected genes (log2-CPM > 1)"
        )

    ggsave(
        file.path(FIG, "groups", "detected_genes_by_group.png"),
        detection_plot, width = 7, height = 6, dpi = 300
    )

    # 11. CLINICAL AND MOLECULAR METADATA

    if ("age" %in% names(metadata)) {

        age_data <- data.frame(
            Age = as.numeric(metadata$age),
            Group = GROUP
        )

        age_plot <- ggplot(age_data, aes(Group, Age, fill = Group)) +
            geom_boxplot() +
            geom_jitter(width = 0.1) +
            theme_bw() +
            labs(title = "Age by DMG subgroup", y = "Age")

        ggsave(
            file.path(FIG, "clinical", "age_by_group.png"),
            age_plot, width = 7, height = 5, dpi = 300
        )
    }

    if ("vaf_h3_K27M" %in% names(metadata)) {

        vaf_data <- data.frame(
            VAF = as.numeric(metadata$vaf_h3_K27M),
            Group = GROUP,
            PC1 = pca$x[, 1]
        )

        vaf_plot <- ggplot(vaf_data, aes(Group, VAF, fill = Group)) +
            geom_boxplot() +
            geom_jitter(width = 0.1) +
            theme_bw() +
            labs(title = "H3 K27M VAF by DMG subgroup", y = "VAF")

        ggsave(
            file.path(FIG, "molecular", "VAF_by_group.png"),
            vaf_plot, width = 7, height = 5, dpi = 300
        )

        vaf_pc1 <- ggplot(vaf_data, aes(VAF, PC1, color = Group)) +
            geom_point(size = 3) +
            geom_smooth(method = "lm", se = FALSE) +
            theme_bw() +
            labs(title = "H3 K27M VAF and PC1", x = "VAF", y = "PC1")

        ggsave(
            file.path(FIG, "molecular", "VAF_vs_PC1.png"),
            vaf_pc1, width = 7, height = 5, dpi = 300
        )
    }

    message("Descriptive statistics completed.")
}

# 12. COMMAND LINE

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript Desc.r DATA_DIR OUTPUT_DIR")
}

descriptive_stats(args[1], args[2])
