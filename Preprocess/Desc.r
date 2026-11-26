#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

descriptive_stats <- function(DATA_DIR, OUTPUT_DIR) {

  # CONFIGURATION
  
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  FIG <- file.path(OUTPUT_DIR, "figures")
  
  for (folder in c(
    "expression",
    "samples",
    "clinical",
    "molecular",
    "groups"
  )) {
    dir.create(
      file.path(FIG, folder),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  LOGCPM_FILE <- file.path(DATA_DIR, "normalized_log2CPM.tsv")
  TPM_FILE <- file.path(DATA_DIR, "merged_samples_SYMBOL_tpm.tsv")
  METADATA_FILE <- file.path(DATA_DIR, "Metadata_for_M2_Biomede.txt")
  
  SAMPLE_ID <- "sample"
  GROUP_COLUMN <- "DMG_subgroup"
  TOP_GENES <- 100
  
  # Variables to explore in PCA.
  # Only variables actually present in metadata are used.
  
  PCA_CATEGORICAL_VARIABLES <- c(
    "Sexe",
    "kit",
    "Localisation",
    "data_omics",
    "Mut_type",
    "Mutation_type",
    "mutation_type",
    "Metastases_at_progression",
    "survival",
    "Survival",
    "progression",
    "Progression",
    "Cancer_progression",
    "cancer_progression",
    "DMG_subtype"
  )
  
  PCA_NUMERIC_VARIABLES <- c(
    "age",
    "vaf_h3_K27M"
  )
  
  # Group colors
  GROUP_COLORS <- c(
    "DMG_H3.3_K27M_TP53wt" = "#377EB8",
    "DMG_H3.3_K27M_TP53mut" = "#E41A1C"
  )
  
  
  # 1. SIGNIFICANCE FUNCTIONS
  
  significance_label <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.0001) return("****")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    
    return("ns")
  }
  
  add_significance <- function(plot, data, value) {
    
    test_data <- data[
      !is.na(data[[value]]) &
        is.finite(data[[value]]) &
        !is.na(data$Group),
      ,
      drop = FALSE
    ]
    
    test_data$Group <- droplevels(factor(test_data$Group))
    
    if (nlevels(test_data$Group) != 2 ||
        any(table(test_data$Group) < 2)) {
      
      warning(
        paste("Cannot test", value,
              ": two groups with at least two observations are required.")
      )
      
      return(list(plot = plot, result = NULL))
    }
    
    test <- wilcox.test(
      test_data[[value]] ~ test_data$Group,
      exact = FALSE
    )
    
    p_value <- test$p.value
    label <- significance_label(p_value)
    
    y_max <- max(test_data[[value]], na.rm = TRUE)
    y_min <- min(test_data[[value]], na.rm = TRUE)
    
    y_range <- y_max - y_min
    
    if (!is.finite(y_range) || y_range == 0) {
      y_range <- max(abs(y_max), 1)
    }
    
    bracket_y <- y_max + 0.12 * y_range
    bracket_height <- 0.035 * y_range
    
    plot <- plot +
      annotate(
        "segment",
        x = 1, xend = 2,
        y = bracket_y,
        yend = bracket_y
      ) +
      annotate(
        "segment",
        x = 1, xend = 1,
        y = bracket_y - bracket_height,
        yend = bracket_y
      ) +
      annotate(
        "segment",
        x = 2, xend = 2,
        y = bracket_y - bracket_height,
        yend = bracket_y
      ) +
      annotate(
        "text",
        x = 1.5,
        y = bracket_y + 0.04 * y_range,
        label = label,
        size = 6
      ) +
      expand_limits(
        y = bracket_y + 0.15 * y_range
      )
    
    result <- data.frame(
      Variable = value,
      Test = "Wilcoxon rank-sum",
      Group_1 = levels(test_data$Group)[1],
      Group_2 = levels(test_data$Group)[2],
      N_Group_1 = sum(
        test_data$Group == levels(test_data$Group)[1]
      ),
      N_Group_2 = sum(
        test_data$Group == levels(test_data$Group)[2]
      ),
      P_value = p_value,
      Significance = label
    )
    
    return(list(
      plot = plot,
      result = result
    ))
  }
  
  significance_results <- list()
  
  
  # 2. LOAD DATA
  
  logCPM <- read.delim(
    LOGCPM_FILE,
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
    na.strings = c("", "NA")
  )
  
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
    match(colnames(logCPM), metadata[[SAMPLE_ID]]),
    ,
    drop = FALSE
  ]
  
  if (anyNA(metadata[[GROUP_COLUMN]])) {
    stop("Missing DMG_subgroup information.")
  }
  
  GROUP <- factor(metadata[[GROUP_COLUMN]])
  
  if (nlevels(GROUP) != 2) {
    stop("Exactly two DMG subgroups are required.")
  }
  
  tpm <- tpm[
    rownames(logCPM),
    colnames(logCPM),
    drop = FALSE
  ]
  
  
  # 3. SAMPLE STATISTICS
  
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
    Detected_Genes_Percent =
      detected / nrow(logCPM) * 100
  )
  
  write.table(
    sample_stats,
    file.path(OUTPUT_DIR, "sample_descriptive_stats.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  group_stats <- aggregate(
    sample_stats[, c(
      "Mean",
      "Median",
      "SD",
      "Detected_Genes",
      "Detected_Genes_Percent"
    )],
    list(Group = GROUP),
    function(x) c(
      mean = mean(x),
      median = median(x),
      sd = sd(x)
    )
  )
  
  write.table(
    group_stats,
    file.path(OUTPUT_DIR, "sample_statistics_by_group.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # 4. GENE STATISTICS
  
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
    order(
      gene_stats$Variance_logCPM,
      decreasing = TRUE
    ),
  ]
  
  write.table(
    gene_stats,
    file.path(OUTPUT_DIR, "gene_descriptive_stats.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    head(gene_stats, TOP_GENES),
    file.path(OUTPUT_DIR, "top_variable_genes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # 5. SAMPLE CORRELATION
  
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
  
  annotation <- data.frame(
    Group = GROUP
  )
  
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
  
  
  # 6. PCA
  
  gene_variances <- apply(logCPM, 1, var)
  selected_genes <- names(
    sort(
      gene_variances,
      decreasing = TRUE
    )
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
    Group = GROUP
  )
  
  write.table(
    pca_coordinates,
    file.path(OUTPUT_DIR, "PCA_coordinates.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    data.frame(
      Component = paste0(
        "PC",
        seq_along(explained)
      ),
      Explained_Variance_Percent = explained
    ),
    file.path(OUTPUT_DIR, "PCA_explained_variance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  pca_plot <- ggplot(
    pca_coordinates,
    aes(
      PC1,
      PC2,
      color = Group
    )
  ) +
    geom_point(size = 3) +
    scale_color_manual(
      values = GROUP_COLORS,
      drop = FALSE
    ) +
    labs(
      title = "Principal Component Analysis",
      x = paste0(
        "PC1 (",
        round(explained[1], 2),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(explained[2], 2),
        "%)"
      ),
      color = "TP53 subgroup"
    ) +
    theme_bw()
  
  ggsave(
    file.path(FIG, "samples", "PCA.png"),
    pca_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # 7. PCA BY CLINICAL AND MOLECULAR METADATA
  
  pca_variables <- intersect(
    c(
      PCA_CATEGORICAL_VARIABLES,
      PCA_NUMERIC_VARIABLES
    ),
    names(metadata)
  )
  
  for (variable in pca_variables) {
    
    values <- as.character(metadata[[variable]])
    
    values[
      values %in% c("", ".", "NA", "NaN")
    ] <- NA
    
    pca_coordinates$Metadata <- values
    
    # Numeric variables: continuous color scale
    
    if (variable %in% PCA_NUMERIC_VARIABLES) {
      
      numeric_values <- suppressWarnings(
        as.numeric(values)
      )
      
      pca_coordinates$Metadata_numeric <- numeric_values
      
      plot <- ggplot(
        pca_coordinates,
        aes(
          PC1,
          PC2,
          color = Metadata_numeric,
          shape = Group
        )
      ) +
        geom_point(
          size = 3,
          alpha = 0.85
        ) +
        scale_color_viridis_c(
          na.value = "grey75"
        ) +
        labs(
          title = paste(
            "PCA colored by",
            variable
          ),
          x = paste0(
            "PC1 (",
            round(explained[1], 2),
            "%)"
          ),
          y = paste0(
            "PC2 (",
            round(explained[2], 2),
            "%)"
          ),
          color = variable,
          shape = "TP53 subgroup"
        ) +
        theme_bw()
      
    } else {
      
      # Categorical variables: discrete color scale
      
      pca_coordinates$Metadata <- factor(
        pca_coordinates$Metadata
      )
      
      plot <- ggplot(
        pca_coordinates,
        aes(
          PC1,
          PC2,
          color = Metadata,
          shape = Group
        )
      ) +
        geom_point(
          size = 3,
          alpha = 0.85
        ) +
        labs(
          title = paste(
            "PCA colored by",
            variable
          ),
          x = paste0(
            "PC1 (",
            round(explained[1], 2),
            "%)"
          ),
          y = paste0(
            "PC2 (",
            round(explained[2], 2),
            "%)"
          ),
          color = variable,
          shape = "TP53 subgroup"
        ) +
        theme_bw()
    }
    
    ggsave(
      file.path(
        FIG,
        "samples",
        paste0("PCA_", variable, ".png")
      ),
      plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  
  # 8. SAMPLE DENDROGRAM WITH TP53 SUBGROUP ANNOTATION
  
  hc <- hclust(
    as.dist(1 - sample_cor),
    method = "average"
  )
  
  dendrogram_groups <- as.character(GROUP)[hc$order]
  
  dendrogram_colors <- unname(
    GROUP_COLORS[dendrogram_groups]
  )
  
  png(
    file.path(FIG, "samples", "sample_dendrogram.png"),
    width = 2200,
    height = 1500,
    res = 200
  )
  
  layout(
    matrix(c(1, 2), nrow = 2),
    heights = c(4, 0.5)
  )
  
  par(mar = c(2, 4, 3, 1))
  
  plot(
    hc,
    main = "Sample hierarchical clustering",
    xlab = "",
    sub = "",
    ylab = "Distance (1 - Pearson correlation)",
    cex = 0.7
  )
  
  legend(
    "topright",
    legend = c("WT", "MUT"),
    fill = c(
      GROUP_COLORS["DMG_H3.3_K27M_TP53wt"],
      GROUP_COLORS["DMG_H3.3_K27M_TP53mut"]
    ),
    border = NA,
    bty = "n"
  )
  
  par(mar = c(0, 4, 1, 1))
  
  plot.new()
  
  plot.window(
    xlim = c(0.5, length(hc$order) + 0.5),
    ylim = c(0, 1)
  )
  
  rect(
    xleft = seq_along(hc$order) - 0.5,
    ybottom = 0,
    xright = seq_along(hc$order) + 0.5,
    ytop = 1,
    col = dendrogram_colors,
    border = NA
  )
  
  dev.off()
  
  
  # 9. EXPRESSION DISTRIBUTION BY DMG SUBGROUP
  
  # One median log2-CPM value per sample.
  # Samples are the independent observations.
  
  expression_sample <- data.frame(
    Sample = colnames(logCPM),
    Group = GROUP,
    Median_log2CPM = apply(
      logCPM,
      2,
      median
    )
  )
  
  expression_plot <- ggplot(
    expression_sample,
    aes(
      Group,
      Median_log2CPM,
      fill = Group
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.12,
      size = 2.5,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = GROUP_COLORS
    ) +
    theme_bw() +
    labs(
      title = "Expression distribution by DMG subgroup",
      x = "DMG subgroup",
      y = "Median log2-CPM per sample"
    )
  
  expression_result <- add_significance(
    expression_plot,
    expression_sample,
    "Median_log2CPM"
  )
  
  expression_plot <- expression_result$plot
  
  if (!is.null(expression_result$result)) {
    significance_results[[length(significance_results) + 1]] <-
      expression_result$result
  }
  
  ggsave(
    file.path(
      FIG,
      "groups",
      "expression_distribution_by_group.png"
    ),
    expression_plot,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  
  # 10. EXPRESSION DENSITY

  expression_df <- data.frame(
    Expression = as.vector(logCPM),
    Sample = rep(
      colnames(logCPM),
      each = nrow(logCPM)
    ),
    Group = rep(
      GROUP,
      each = nrow(logCPM)
    )
  )
  
  density_plot <- ggplot(
    expression_df,
    aes(
      Expression,
      group = Sample,
      color = Group
    )
  ) +
    geom_density(alpha = 0.15) +
    scale_color_manual(
      values = GROUP_COLORS
    ) +
    theme_bw() +
    labs(
      title = "Expression density across samples",
      x = "log2-CPM",
      y = "Density",
      color = "TP53 subgroup"
    )
  
  ggsave(
    file.path(
      FIG,
      "expression",
      "expression_density.png"
    ),
    density_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # 11. MEAN-VARIANCE PLOT
  
  mean_variance <- ggplot(
    gene_stats,
    aes(
      Mean_logCPM,
      Variance_logCPM
    )
  ) +
    geom_point(
      alpha = 0.3,
      size = 0.7
    ) +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = "Mean-variance relationship",
      x = "Mean log2-CPM",
      y = "Variance of log2-CPM"
    )
  
  ggsave(
    file.path(
      FIG,
      "expression",
      "mean_variance.png"
    ),
    mean_variance,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # 12. DETECTED GENES BY GROUP
  
  detection_plot <- ggplot(
    sample_stats,
    aes(
      Group,
      Detected_Genes,
      fill = Group
    )
  ) +
    geom_boxplot() +
    geom_jitter(
      width = 0.1,
      size = 2
    ) +
    scale_fill_manual(
      values = GROUP_COLORS
    ) +
    theme_bw() +
    labs(
      title = "Detected genes by DMG subgroup",
      x = "DMG subgroup",
      y = "Number of detected genes (log2-CPM > 1)"
    )
  
  detection_result <- add_significance(
    detection_plot,
    sample_stats,
    "Detected_Genes"
  )
  
  detection_plot <- detection_result$plot
  
  if (!is.null(detection_result$result)) {
    significance_results[[length(significance_results) + 1]] <-
      detection_result$result
  }
  
  ggsave(
    file.path(
      FIG,
      "groups",
      "detected_genes_by_group.png"
    ),
    detection_plot,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  
  # 13. AGE BY DMG SUBGROUP
  
  if ("age" %in% names(metadata)) {
    
    age_data <- data.frame(
      Age = suppressWarnings(
        as.numeric(metadata$age)
      ),
      Group = GROUP
    )
    
    age_plot <- ggplot(
      age_data,
      aes(
        Group,
        Age,
        fill = Group
      )
    ) +
      geom_boxplot() +
      geom_jitter(
        width = 0.1,
        size = 2
      ) +
      scale_fill_manual(
        values = GROUP_COLORS
      ) +
      theme_bw() +
      labs(
        title = "Age by DMG subgroup",
        x = "DMG subgroup",
        y = "Age"
      )
    
    age_result <- add_significance(
      age_plot,
      age_data,
      "Age"
    )
    
    age_plot <- age_result$plot
    
    if (!is.null(age_result$result)) {
      significance_results[[length(significance_results) + 1]] <-
        age_result$result
    }
    
    ggsave(
      file.path(
        FIG,
        "clinical",
        "age_by_group.png"
      ),
      age_plot,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  
  # 14. H3 K27M VAF
  
  if ("vaf_h3_K27M" %in% names(metadata)) {
    
    
    # VAF DATA
    
    vaf_data <- data.frame(
      VAF = suppressWarnings(
        as.numeric(as.character(metadata$vaf_h3_K27M))
      ),
      Group = GROUP,
      PC1 = pca$x[, 1],
      Name = as.character(metadata$sample)
    )
    
    # Remove non-finite VAF values from outlier detection
    vaf_data$Outlier <- FALSE
    
    
    # OUTLIER DETECTION
    
    vaf_data$Outlier <- FALSE
    
    message("\n========== VAF OUTLIER DIAGNOSTIC ==========")
    
    for (grp in levels(GROUP)) {
      
      idx <- which(
        vaf_data$Group == grp &
          is.finite(vaf_data$VAF)
      )
      
      if (length(idx) == 0) next
      
      values <- vaf_data$VAF[idx]
      
      box_stats <- boxplot.stats(values)
      
      lower_whisker <- box_stats$stats[1]
      upper_whisker <- box_stats$stats[5]
      
      lower_outliers <- values[values < lower_whisker]
      upper_outliers <- values[values > upper_whisker]
      
      # Identify every value outside the box
      lower_box <- box_stats$stats[2]
      upper_box <- box_stats$stats[4]
      
      vaf_data$Outlier[idx] <-
        vaf_data$VAF[idx] < lower_box |
        vaf_data$VAF[idx] > upper_box
      
      message("\nGROUP: ", grp)
      
      message(
        "Lower whisker: ", round(lower_whisker, 4),
        "\nUpper whisker: ", round(upper_whisker, 4),
        "\nLower outliers: ", length(lower_outliers),
        "\nUpper outliers: ", length(upper_outliers),
        "\nTotal outliers: ", length(box_stats$out)
      )
    }
    
    
    # OUTLIER TABLE
    
    
    vaf_outliers <- vaf_data[
      vaf_data$Outlier &
        !is.na(vaf_data$Name) &
        nzchar(vaf_data$Name),
      ,
      drop = FALSE
    ]
    
    message(
      "\nTotal VAF outliers detected: ",
      nrow(vaf_outliers)
    )
    
    if (nrow(vaf_outliers) > 0) {
      print(
        vaf_outliers[, c("Name", "Group", "VAF")]
      )
    }
    
    
    # POINT POSITIONS
    
    
    set.seed(123)
    
    vaf_data$Group_x <- as.numeric(vaf_data$Group)
    
    vaf_data$Point_x <- vaf_data$Group_x +
      runif(nrow(vaf_data), -0.08, 0.08)
    
    vaf_outliers$Group_x <- as.numeric(vaf_outliers$Group)
    
    vaf_outliers$Point_x <- vaf_data$Point_x[
      match(vaf_outliers$Name, vaf_data$Name)
    ]
    
    
    # VAF BOXPLOT
    
    
    vaf_plot <- ggplot() +
      
      geom_boxplot(
        data = vaf_data,
        aes(
          x = Group_x,
          y = VAF,
          group = Group,
          fill = Group
        ),
        outlier.shape = NA
      ) +
      
      geom_point(
        data = vaf_data,
        aes(
          x = Point_x,
          y = VAF
        ),
        color = "black",
        size = 2
      ) +
      
      # Label only detected outliers
      ggrepel::geom_text_repel(
        data = vaf_outliers,
        aes(
          x = Point_x,
          y = VAF,
          label = Name
        ),
        direction = "y",
        nudge_x = 0.15,
        hjust = 0,
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        min.segment.length = 0,
        seed = 123
      ) +
      
      scale_x_continuous(
        breaks = seq_along(levels(GROUP)),
        labels = levels(GROUP)
      ) +
      
      scale_fill_manual(
        values = GROUP_COLORS
      ) +
      
      coord_cartesian(
        clip = "off"
      ) +
      
      theme_bw() +
      
      theme(
        plot.margin = margin(
          5.5, 60, 5.5, 5.5
        )
      ) +
      
      labs(
        title = "H3 K27M VAF by DMG subgroup",
        x = "DMG subgroup",
        y = "VAF"
      )
    
    
    # SIGNIFICANCE TEST
    
    
    vaf_result <- add_significance(
      vaf_plot,
      vaf_data,
      "VAF"
    )
    
    vaf_plot <- vaf_result$plot
    
    if (!is.null(vaf_result$result)) {
      significance_results[[length(significance_results) + 1]] <-
        vaf_result$result
    }
    
    
    # SAVE VAF BOXPLOT
    
    
    ggsave(
      file.path(
        FIG,
        "molecular",
        "VAF_by_group.png"
      ),
      vaf_plot,
      width = 10,
      height = 5,
      dpi = 300
    )
    
    
    # VAF VS PC1
    
    
    vaf_pc1 <- ggplot(
      vaf_data,
      aes(
        VAF,
        PC1,
        color = Group
      )
    ) +
      
      geom_point(
        size = 3
      ) +
      
      geom_smooth(
        method = "lm",
        se = FALSE
      ) +
      
      scale_color_manual(
        values = GROUP_COLORS
      ) +
      
      theme_bw() +
      
      labs(
        title = "H3 K27M VAF and PC1",
        x = "VAF",
        y = "PC1",
        color = "TP53 subgroup"
      )
    
    ggsave(
      file.path(
        FIG,
        "molecular",
        "VAF_vs_PC1.png"
      ),
      vaf_pc1,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  

  # 15. EMT GENE EXPRESSION

  EMT_GENES <- c(
    "SNAI1",
    "TWIST1",
    "FN1",
    "CDH2",
    "CDH11"
  )
  
  # Identify genes available after expression filtering
  
  available_genes <- intersect(
    EMT_GENES,
    rownames(logCPM)
  )
  
  missing_genes <- setdiff(
    EMT_GENES,
    available_genes
  )
  
  if (length(missing_genes) > 0) {
    warning(
      paste(
        "EMT genes not found in expression matrix:",
        paste(missing_genes, collapse = ", ")
      )
    )
  }
  
  if (length(available_genes) > 0) {
    
    # Prepare expression data
    
    emt_data <- do.call(
      rbind,
      lapply(available_genes, function(gene) {
        
        data.frame(
          Gene = gene,
          Sample = colnames(logCPM),
          Expression = as.numeric(logCPM[gene, ]),
          Group = GROUP
        )
      })
    )
    
    emt_data$Gene <- factor(
      emt_data$Gene,
      levels = available_genes
    )
    
    # Wilcoxon test for each gene
    
    emt_results <- lapply(
      available_genes,
      function(gene) {
        
        gene_data <- emt_data[
          emt_data$Gene == gene &
            is.finite(emt_data$Expression),
          ,
          drop = FALSE
        ]
        
        gene_data$Group <- droplevels(
          factor(gene_data$Group)
        )
        
        if (
          nlevels(gene_data$Group) != 2 ||
          any(table(gene_data$Group) < 2)
        ) {
          return(data.frame(
            Gene = gene,
            N_WT = NA,
            N_MUT = NA,
            Median_WT = NA,
            Median_MUT = NA,
            P_value = NA
          ))
        }
        
        test <- wilcox.test(
          Expression ~ Group,
          data = gene_data,
          exact = FALSE
        )
        
        data.frame(
          Gene = gene,
          N_WT = sum(
            gene_data$Group ==
              "DMG_H3.3_K27M_TP53wt"
          ),
          N_MUT = sum(
            gene_data$Group ==
              "DMG_H3.3_K27M_TP53mut"
          ),
          Median_WT = median(
            gene_data$Expression[
              gene_data$Group ==
                "DMG_H3.3_K27M_TP53wt"
            ]
          ),
          Median_MUT = median(
            gene_data$Expression[
              gene_data$Group ==
                "DMG_H3.3_K27M_TP53mut"
            ]
          ),
          P_value = test$p.value
        )
      }
    )
    
    emt_results <- do.call(
      rbind,
      emt_results
    )
    
    # Multiple testing correction across EMT genes
    
    emt_results$Adjusted_P_value <- p.adjust(
      emt_results$P_value,
      method = "BH"
    )
    
    emt_results$Significance <- vapply(
      emt_results$Adjusted_P_value,
      significance_label,
      character(1)
    )
    
    write.table(
      emt_results,
      file.path(
        OUTPUT_DIR,
        "EMT_gene_statistics.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # Prepare significance annotations
    
    annotation_data <- lapply(
      available_genes,
      function(gene) {
        
        gene_values <- emt_data$Expression[
          emt_data$Gene == gene
        ]
        
        gene_values <- gene_values[
          is.finite(gene_values)
        ]
        
        y_max <- max(gene_values)
        y_min <- min(gene_values)
        y_range <- y_max - y_min
        
        if (y_range == 0) {
          y_range <- max(abs(y_max), 1)
        }
        
        result <- emt_results[
          emt_results$Gene == gene,
          ,
          drop = FALSE
        ]
        
        data.frame(
          Gene = gene,
          y = y_max + 0.15 * y_range,
          label = result$Significance
        )
      }
    )
    
    annotation_data <- do.call(
      rbind,
      annotation_data
    )
    
    annotation_data$Gene <- factor(
      annotation_data$Gene,
      levels = available_genes
    )
    
    # Plot: one panel per gene, WT vs MUT
    
    emt_plot <- ggplot(
      emt_data,
      aes(
        x = Group,
        y = Expression,
        fill = Group
      )
    ) +
      geom_boxplot(
        width = 0.55,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.12,
        size = 1.8,
        color = "black",
        alpha = 0.75,
        show.legend = FALSE
      ) +
      facet_wrap(
        ~ Gene,
        scales = "free_y",
        ncol = 3
      ) +
      geom_segment(
        data = annotation_data,
        aes(
          x = 1,
          xend = 2,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = annotation_data,
        aes(
          x = 1,
          xend = 1,
          y = y - 0.04 * abs(y),
          yend = y
        ),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = annotation_data,
        aes(
          x = 2,
          xend = 2,
          y = y - 0.04 * abs(y),
          yend = y
        ),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = annotation_data,
        aes(
          x = 1.5,
          y = y + 0.04 * abs(y),
          label = label
        ),
        inherit.aes = FALSE,
        size = 5
      ) +
      scale_fill_manual(
        values = GROUP_COLORS
      ) +
      scale_color_manual(
        values = GROUP_COLORS
      ) +
      scale_y_continuous(
        expand = expansion(
          mult = c(0.05, 0.25)
        )
      ) +
      labs(
        title = "EMT gene expression in H3K27M gliomas",
        subtitle = "TP53 WT vs MUT",
        x = "TP53 subgroup",
        y = "log2-CPM",
        fill = "TP53 subgroup"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(
          face = "bold",
          size = 12
        ),
        axis.text.x = element_text(
          angle = 20,
          hjust = 1
        )
      )
    
    ggsave(
      file.path(
        FIG,
        "molecular",
        "EMT_gene_expression.png"
      ),
      emt_plot,
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  
  # 16. DIFFERENTIAL EXPRESSION AND VOLCANO PLOT

  # Define TP53 groups explicitly
  
  GROUP_SHORT <- factor(
    ifelse(
      grepl("mut$", as.character(GROUP), ignore.case = TRUE),
      "MUT",
      "WT"
    ),
    levels = c("WT", "MUT")
  )
  
  # Design matrix
  
  design <- model.matrix(~ 0 + GROUP_SHORT)
  
  colnames(design) <- levels(GROUP_SHORT)
  
  # Fit linear model on normalized log2-CPM
  
  fit <- limma::lmFit(
    logCPM,
    design
  )
  
  # Compare MUT against WT
  
  contrast <- limma::makeContrasts(
    MUT - WT,
    levels = design
  )
  
  fit_contrast <- limma::contrasts.fit(
    fit,
    contrast
  )
  
  # Empirical Bayes moderation
  
  fit_contrast <- limma::eBayes(
    fit_contrast,
    trend = TRUE,
    robust = TRUE
  )
  
  # Extract results for all genes
  
  volcano_results <- limma::topTable(
    fit_contrast,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  volcano_results$Gene <- rownames(volcano_results)
  
  # Save complete differential expression table
  
  write.table(
    volcano_results,
    file.path(
      OUTPUT_DIR,
      "volcano_differential_expression.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Significance thresholds
  
  LOGFC_THRESHOLD <- 1
  PADJ_THRESHOLD <- 0.05
  
  volcano_results$Significance <- "Not significant"
  
  volcano_results$Significance[
    volcano_results$logFC >= LOGFC_THRESHOLD &
      volcano_results$adj.P.Val < PADJ_THRESHOLD
  ] <- "Upregulated in MUT"
  
  volcano_results$Significance[
    volcano_results$logFC <= -LOGFC_THRESHOLD &
      volcano_results$adj.P.Val < PADJ_THRESHOLD
  ] <- "Downregulated in MUT"
  
  volcano_results$Significance <- factor(
    volcano_results$Significance,
    levels = c(
      "Downregulated in MUT",
      "Not significant",
      "Upregulated in MUT"
    )
  )
  
  # Avoid infinite values when adjusted p-values are zero
  
  volcano_results$NegLog10_Padj <- -log10(
    pmax(volcano_results$adj.P.Val, 1e-300)
  )
  
  # Label the 10 most significant genes
  
  significant_genes <- volcano_results[
    volcano_results$adj.P.Val < PADJ_THRESHOLD &
      abs(volcano_results$logFC) >= LOGFC_THRESHOLD,
    ,
    drop = FALSE
  ]
  
  significant_genes <- significant_genes[
    order(significant_genes$adj.P.Val),
    ,
    drop = FALSE
  ]
  
  labels_to_display <- head(
    significant_genes,
    10
  )
  
  # Volcano plot
  
  volcano_plot <- ggplot(
    volcano_results,
    aes(
      x = logFC,
      y = NegLog10_Padj,
      color = Significance
    )
  ) +
    geom_point(
      alpha = 0.7,
      size = 1.5
    ) +
    geom_vline(
      xintercept = c(
        -LOGFC_THRESHOLD,
        LOGFC_THRESHOLD
      ),
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(PADJ_THRESHOLD),
      linetype = "dashed",
      color = "black"
    ) +
    geom_text(
      data = labels_to_display,
      aes(label = Gene),
      inherit.aes = TRUE,
      size = 3,
      vjust = -0.7,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Downregulated in MUT" = "#377EB8",
        "Not significant" = "grey75",
        "Upregulated in MUT" = "#E41A1C"
      )
    ) +
    theme_bw() +
    labs(
      title = "Differential expression: TP53 MUT vs WT",
      subtitle = paste(
        "log2FC threshold:",
        LOGFC_THRESHOLD,
        "| Adjusted p-value:",
        PADJ_THRESHOLD
      ),
      x = "log2 fold change (MUT - WT)",
      y = "-log10 adjusted p-value (BH)",
      color = "Expression status"
    )
  
  ggsave(
    file.path(
      FIG,
      "groups",
      "volcano_plot.png"
    ),
    volcano_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  # Export significant genes separately
  
  write.table(
    significant_genes,
    file.path(
      OUTPUT_DIR,
      "significant_DE_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message(
    "Volcano plot completed. Significant genes: ",
    nrow(significant_genes)
  )
  
  
  # 17. SIGNIFICANCE RESULTS
  
  if (length(significance_results) > 0) {
    
    significance_table <- do.call(
      rbind,
      significance_results
    )
    
    significance_table$Adjusted_P_value <- p.adjust(
      significance_table$P_value,
      method = "BH"
    )
    
    write.table(
      significance_table,
      file.path(
        OUTPUT_DIR,
        "group_significance_tests.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  message("Descriptive statistics completed.")
}


# COMMAND LINE

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript Desc.r DATA_DIR OUTPUT_DIR")
}

descriptive_stats(
  args[1],
  args[2]
)
