#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(scales)
})

if (file.exists("R/_shared.R")) source("R/_shared.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

dir.create("results/supplementary_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_package/supplementary_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.45, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_pub <- function(p, stem, width = 8, height = 6, dpi = 450) {
  for (d in c("results/supplementary_figures", "paper_package/supplementary_figures")) {
    ggsave(file.path(d, paste0(stem, ".pdf")), p,
           width = width, height = height, units = "in",
           dpi = dpi, limitsize = FALSE, bg = "white")
    ggsave(file.path(d, paste0(stem, ".png")), p,
           width = width, height = height, units = "in",
           dpi = dpi, limitsize = FALSE, bg = "white")
  }
  message("Wrote: ", stem)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_col <- function(df, patterns, required = TRUE, label = "column") {
  hits <- unique(unlist(lapply(patterns, function(p) {
    grep(p, colnames(df), value = TRUE, ignore.case = TRUE)
  })))
  if (length(hits) == 0) {
    if (required) {
      stop("Could not find ", label, ". Columns: ", paste(colnames(df), collapse = ", "))
    }
    return(NA_character_)
  }
  hits[1]
}

combine_plots <- function(..., ncol = 1, nrow = NULL) {
  plots <- list(...)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(Reduce(`+`, plots) + patchwork::plot_layout(ncol = ncol, nrow = nrow))
  }
  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(plotlist = plots, ncol = ncol, nrow = nrow, align = "hv"))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow))
  }
  stop("Install one of: patchwork, cowplot, or gridExtra.")
}

load_curated <- function(acc) {
  f <- file.path("data/processed", paste0(acc, "_curated.rds"))
  if (!file.exists(f)) stop("Missing curated object: ", f)
  readRDS(f)
}

# ---------------------------
# Supplementary Figure S1
# MDS quality control
# ---------------------------

make_s1 <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package limma is required.")
  
  mds_df <- function(acc, color_by) {
    obj <- load_curated(acc)
    expr <- obj$expr
    meta <- obj$metadata |> mutate(sample_id = as.character(sample_id))
    
    md <- limma::plotMDS(expr, plot = FALSE)
    
    df <- tibble(
      sample_id = colnames(expr),
      Dim1 = md$x,
      Dim2 = md$y
    ) |>
      left_join(meta, by = "sample_id")
    
    if (!(color_by %in% colnames(df))) {
      geno_col <- find_col(df, c("IL28B", "rs8099917", "genotype"), required = FALSE)
      df[[color_by]] <- if (!is.na(geno_col)) df[[geno_col]] else NA_character_
    }
    
    df |>
      mutate(
        dataset = acc,
        color_by = color_by,
        color_value = as.character(.data[[color_by]])
      )
  }
  
  plot_mds <- function(df, title, legend_title) {
    ggplot(df, aes(Dim1, Dim2, fill = color_value)) +
      geom_point(shape = 21, size = 2.8, color = "black", stroke = 0.25, alpha = 0.9) +
      labs(title = title, x = "MDS dimension 1", y = "MDS dimension 2", fill = legend_title) +
      theme_pub(10)
  }
  
  d1 <- mds_df("GSE121248", "tissue")
  d2 <- mds_df("GSE41804", "tissue")
  d3 <- mds_df("GSE41804", "IL28B_genotype")
  
  p <- combine_plots(
    plot_mds(d1, "GSE121248: tissue structure", "Tissue"),
    plot_mds(d2, "GSE41804: tissue structure", "Tissue"),
    plot_mds(d3, "GSE41804: IL28B genotype", "IL28B genotype"),
    ncol = 1
  )
  
  save_pub(p, "Supplementary_Figure_S1_MDS_QC", width = 7.2, height = 10.5)
  
  write_tsv(bind_rows(d1, d2, d3), "results/tables/suppfig_s1_mds_plot_data.tsv")
}

# ---------------------------
# Supplementary Figure S2
# HBV Hallmark volcano + heatmap
# ---------------------------

make_s2 <- function() {
  limma_file <- first_existing(c(
    "results/tables/GSE121248_hallmark_limma.tsv",
    "results/tables/GSE121248_hallmark_results.tsv",
    "results/tables/hallmark_GSE121248_limma.tsv"
  ))
  
  if (is.na(limma_file)) stop("Could not find GSE121248 Hallmark limma table.")
  
  hall <- read_tsv(limma_file, show_col_types = FALSE)
  
  set_col <- find_col(hall, c("Hallmark", "hallmark", "gs_name", "pathway", "gene_set"))
  logfc_col <- find_col(hall, c("^logFC$", "log_fc", "estimate"))
  p_col <- find_col(hall, c("^P.Value$", "^p_value$", "^p$", "pval"))
  fdr_col <- find_col(hall, c("^FDR$", "adj.P.Val", "adj_p", "padj"), required = FALSE)
  
  hall2 <- hall |>
    transmute(
      pathway_raw = as.character(.data[[set_col]]),
      pathway = str_remove(pathway_raw, "^HALLMARK_"),
      logFC = as.numeric(.data[[logfc_col]]),
      p_value = as.numeric(.data[[p_col]]),
      FDR = if (!is.na(fdr_col)) as.numeric(.data[[fdr_col]]) else p.adjust(p_value, "BH")
    ) |>
    mutate(
      neglog10p = -log10(pmax(p_value, .Machine$double.xmin)),
      direction = case_when(
        FDR < 0.05 & logFC > 0 ~ "Activated in tumor",
        FDR < 0.05 & logFC < 0 ~ "Suppressed in tumor",
        TRUE ~ "Not FDR-significant"
      )
    )
  
  volcano <- ggplot(hall2, aes(logFC, neglog10p, fill = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey40") +
    geom_point(shape = 21, size = 2.8, color = "black", stroke = 0.22, alpha = 0.92) +
    scale_fill_manual(values = c(
      "Activated in tumor" = "#D55E00",
      "Suppressed in tumor" = "#0072B2",
      "Not FDR-significant" = "grey75"
    )) +
    labs(
      title = "GSE121248 Hallmark differential activity",
      x = "Tumor vs adjacent/non-tumor logFC",
      y = expression(-log[10](p)),
      fill = NULL
    ) +
    theme_pub(10) +
    theme(legend.position = "bottom")
  
  score_file <- first_existing(c(
    "data/processed/GSE121248_hallmark_scores.rds",
    "results/gsva/GSE121248_hallmark_scores.rds"
  ))
  
  if (is.na(score_file)) {
    save_pub(volcano, "Supplementary_Figure_S2_HBV_Hallmark_volcano", width = 7.5, height = 5.3)
    write_tsv(hall2, "results/tables/suppfig_s2_hbv_hallmark_volcano_data.tsv")
    return(invisible(NULL))
  }
  
  score_obj <- readRDS(score_file)
  scores <- if (is.list(score_obj) && "scores" %in% names(score_obj)) score_obj$scores else score_obj
  scores <- as.matrix(scores)
  
  meta <- load_curated("GSE121248")$metadata |>
    mutate(sample_id = as.character(sample_id))
  
  top_paths <- bind_rows(
    hall2 |> filter(logFC > 0) |> arrange(FDR, desc(logFC)) |> slice_head(n = 12),
    hall2 |> filter(logFC < 0) |> arrange(FDR, logFC) |> slice_head(n = 12)
  ) |>
    pull(pathway) |>
    unique()
  
  row_map <- tibble(
    raw = rownames(scores),
    clean = str_remove(rownames(scores), "^HALLMARK_")
  ) |>
    filter(clean %in% top_paths)
  
  mat <- scores[row_map$raw, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  meta2 <- tibble(sample_id = colnames(mat_z)) |>
    left_join(meta, by = "sample_id") |>
    mutate(tissue = factor(tissue, levels = c("non_tumor", "tumor")))
  
  sample_order <- meta2 |>
    arrange(tissue, sample_id) |>
    pull(sample_id)
  
  mat_z <- mat_z[, sample_order, drop = FALSE]
  
  pathway_order <- hall2 |>
    filter(pathway %in% row_map$clean) |>
    arrange(desc(logFC)) |>
    pull(pathway)
  
  raw_order <- row_map$raw[match(pathway_order, row_map$clean)]
  mat_z <- mat_z[raw_order, , drop = FALSE]
  rownames(mat_z) <- str_replace_all(pathway_order, "_", " ")
  
  heat_df <- as.data.frame(mat_z) |>
    rownames_to_column("pathway") |>
    pivot_longer(-pathway, names_to = "sample_id", values_to = "z_score") |>
    left_join(meta2 |> select(sample_id, tissue), by = "sample_id") |>
    mutate(
      sample_id = factor(sample_id, levels = sample_order),
      pathway = factor(pathway, levels = rev(rownames(mat_z)))
    )
  
  heatmap_plot <- ggplot(heat_df, aes(sample_id, pathway, fill = z_score)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Row z-score",
      oob = squish
    ) +
    facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
    labs(title = "GSE121248 top Hallmark activity patterns", x = "Samples", y = NULL) +
    theme_pub(8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )
  
  p <- combine_plots(volcano, heatmap_plot, ncol = 1)
  
  save_pub(p, "Supplementary_Figure_S2_HBV_Hallmark_details", width = 8.5, height = 10.5)
  
  write_tsv(hall2, "results/tables/suppfig_s2_hbv_hallmark_volcano_data.tsv")
  write_tsv(heat_df, "results/tables/suppfig_s2_hbv_hallmark_heatmap_data.tsv")
}

# ---------------------------
# Supplementary Figure S3
# Module-size robustness
# ---------------------------

make_s3 <- function() {
  table_file <- first_existing(c(
    "results/tables/module_size_robustness.tsv",
    "results/tables/module_size_robustness_summary.tsv",
    "results/tables/hccstate_module_size_robustness.tsv",
    "results/tables/validation_module_size_robustness.tsv",
    "results/tables/module_size_sensitivity.tsv"
  ))
  
  if (is.na(table_file)) {
    hits <- list.files(
      "results/tables",
      pattern = "module.*size|size.*robust|module.*sensitivity",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(hits) > 0) table_file <- hits[1]
  }
  
  if (is.na(table_file)) {
    stop("Could not find module-size robustness table. Run the module-size sensitivity pipeline first.")
  }
  
  dat <- read_tsv(table_file, show_col_types = FALSE)
  
  dataset_col <- find_col(dat, c("^dataset$", "^cohort$", "^gse$"))
  topn_col <- find_col(dat, c("^topN$", "top_n", "module_size", "n_genes", "size"))
  auc_col <- find_col(dat, c("^AUC$", "auc"), required = FALSE)
  delta_col <- find_col(dat, c("delta", "difference", "tumor_minus", "estimate"), required = FALSE)
  score_col <- find_col(dat, c("^score$", "module", "metric"), required = FALSE)
  
  plot_dat <- dat |>
    mutate(
      dataset = as.character(.data[[dataset_col]]),
      topN = as.numeric(.data[[topn_col]]),
      score = if (!is.na(score_col)) as.character(.data[[score_col]]) else "HCCStateScore",
      AUC = if (!is.na(auc_col)) as.numeric(.data[[auc_col]]) else NA_real_,
      delta = if (!is.na(delta_col)) as.numeric(.data[[delta_col]]) else NA_real_
    ) |>
    filter(str_detect(score, regex("HCCState|HCC", TRUE)) | score == "HCCStateScore")
  
  plots <- list()
  
  if (!all(is.na(plot_dat$AUC))) {
    plots[[length(plots) + 1]] <- plot_dat |>
      filter(!is.na(AUC)) |>
      ggplot(aes(topN, AUC, group = dataset, color = dataset)) +
      geom_line(linewidth = 0.55) +
      geom_point(size = 2.0) +
      scale_x_continuous(breaks = sort(unique(plot_dat$topN))) +
      coord_cartesian(
        ylim = c(max(0.5, min(plot_dat$AUC, na.rm = TRUE) - 0.03), 1.01),
        clip = "off"
      ) +
      labs(
        title = "HCCStateScore discrimination across module sizes",
        x = "Genes per module",
        y = "AUC",
        color = "Dataset"
      ) +
      theme_pub(10)
  }
  
  if (!all(is.na(plot_dat$delta))) {
    plots[[length(plots) + 1]] <- plot_dat |>
      filter(!is.na(delta)) |>
      ggplot(aes(topN, delta, group = dataset, color = dataset)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.3) +
      geom_line(linewidth = 0.55) +
      geom_point(size = 2.0) +
      scale_x_continuous(breaks = sort(unique(plot_dat$topN))) +
      labs(
        title = "HCCStateScore tumor/non-tumor delta across module sizes",
        x = "Genes per module",
        y = "Tumor - non-tumor delta",
        color = "Dataset"
      ) +
      theme_pub(10)
  }
  
  if (length(plots) == 0) stop("No usable AUC or delta columns found.")
  
  p <- if (length(plots) == 1) plots[[1]] else combine_plots(plots[[1]], plots[[2]], ncol = 1)
  
  save_pub(
    p,
    "Supplementary_Figure_S3_module_size_robustness",
    width = 9.5,
    height = ifelse(length(plots) == 1, 5.8, 9.5)
  )
  
  write_tsv(plot_dat, "results/tables/suppfig_s3_module_size_robustness_plot_data.tsv")
}

# ---------------------------
# Supplementary Figure S4
# HBV injury top-N sensitivity
# ---------------------------

make_s4 <- function() {
  table_file <- first_existing(c(
    "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv",
    "results/tables/gse121248_hbv_injury_topN_extended_regression_tissue.tsv"
  ))
  
  if (is.na(table_file)) stop("Could not find HBV injury top-N summary. Run script 09 first.")
  
  dat <- read_tsv(table_file, show_col_types = FALSE)
  
  set_col <- find_col(dat, c("^injury_set$", "injury"))
  model_col <- find_col(dat, c("^model$"))
  estimate_col <- find_col(dat, c("tumor_coefficient", "^estimate$"))
  low_col <- find_col(dat, c("^ci_low$", "^conf.low$", "conf_low"))
  high_col <- find_col(dat, c("^ci_high$", "^conf.high$", "conf_high"))
  p_col <- find_col(dat, c("^p_value$", "^p.value$", "^p$"))
  n_col <- find_col(dat, c("n_input_genes", "n_genes"), required = FALSE)
  
  model_labels <- c(
    "unadjusted" = "Unadjusted",
    "proliferation_adjusted" = "E2F/G2M adjusted",
    "proliferation_cibersortx_pc_adjusted" = "E2F/G2M + CIBERSORTx PCs",
    "proliferation_selected_fraction_adjusted" = "E2F/G2M + selected fractions"
  )
  
  plot_dat <- dat |>
    transmute(
      injury_set = as.character(.data[[set_col]]),
      n_input_genes = if (!is.na(n_col)) as.numeric(.data[[n_col]]) else NA_real_,
      model = as.character(.data[[model_col]]),
      tumor_coefficient = as.numeric(.data[[estimate_col]]),
      ci_low = as.numeric(.data[[low_col]]),
      ci_high = as.numeric(.data[[high_col]]),
      p_value = as.numeric(.data[[p_col]])
    ) |>
    mutate(
      set_size = case_when(
        str_detect(injury_set, "200$|TOP_200\\b") ~ 200,
        str_detect(injury_set, "500$|TOP_500\\b") ~ 500,
        str_detect(injury_set, "1000$|TOP_1000\\b") ~ 1000,
        str_detect(injury_set, "2000$|TOP_2000\\b") ~ 2000,
        str_detect(injury_set, "5000$|TOP_5000\\b") ~ 5000,
        str_detect(injury_set, "7792|EXTENDED") ~ 7792,
        !is.na(n_input_genes) ~ n_input_genes,
        TRUE ~ NA_real_
      ),
      set_label = case_when(
        set_size == 7792 ~ "Extended 7792",
        !is.na(set_size) ~ paste0("Top ", set_size),
        TRUE ~ injury_set
      ),
      set_label = factor(
        set_label,
        levels = c("Top 200", "Top 500", "Top 1000", "Top 2000", "Top 5000", "Extended 7792")
      ),
      model_label = recode(model, !!!model_labels),
      model_label = factor(model_label, levels = unname(model_labels)),
      significant = p_value < 0.05
    ) |>
    filter(!is.na(set_label), !is.na(tumor_coefficient))
  
  p <- ggplot(
    plot_dat,
    aes(set_label, tumor_coefficient, ymin = ci_low, ymax = ci_high,
        color = model_label, group = model_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_errorbar(position = position_dodge(width = 0.58), width = 0.20, linewidth = 0.45) +
    geom_point(aes(shape = significant), position = position_dodge(width = 0.58), size = 2.4) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), name = "p < 0.05") +
    labs(
      title = "HBV injury top-N sensitivity analysis",
      x = "HBV injury gene-set definition",
      y = "Tumor coefficient",
      color = "Model"
    ) +
    theme_pub(10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right") +
    coord_cartesian(clip = "off")
  
  save_pub(p, "Supplementary_Figure_S4_HBV_INJURY_topN_sensitivity", width = 10.5, height = 6.1)
  
  write_tsv(plot_dat, "results/tables/suppfig_s4_hbv_injury_topn_plot_data.tsv")
}

# ---------------------------
# Supplementary Figure S5
# CIBERSORTx immune composition
# ---------------------------

make_s5 <- function() {
  input_file <- first_existing(c(
    "results/tables/gse121248_hbv_injury_topN_extended_regression_input.tsv",
    "results/tables/gse121248_cibersortx_regression_input.tsv"
  ))
  
  if (is.na(input_file)) stop("Could not find CIBERSORTx regression input. Run script 09 first.")
  
  dat <- read_tsv(input_file, show_col_types = FALSE) |>
    distinct()
  
  if ("injury_set" %in% colnames(dat)) {
    dat <- dat |>
      filter(str_detect(injury_set, "TOP_2000|EXTENDED")) |>
      group_by(sample_id) |>
      slice(1) |>
      ungroup()
  }
  
  pc_cols <- grep("^CIBERSORTx_PC[0-9]+$", colnames(dat), value = TRUE)
  
  exclude <- c(
    "HBV_INJURY", "E2F", "G2M", "ProlifHubScore", "HepLossScore",
    "HCCStateScore", pc_cols, "n_input_genes", "n_overlap_genes"
  )
  
  frac_cols <- setdiff(colnames(dat), c("sample_id", "tissue", "injury_set", "analysis_role", exclude))
  frac_cols <- frac_cols[vapply(dat[frac_cols], is.numeric, logical(1))]
  frac_cols <- frac_cols[!str_detect(frac_cols, regex("score|coefficient|p.value|p_value|ci_|conf|percent|rank", TRUE))]
  
  plots <- list()
  
  if (length(pc_cols) >= 2) {
    plots[[length(plots) + 1]] <- dat |>
      filter(is.finite(.data[[pc_cols[1]]]), is.finite(.data[[pc_cols[2]]])) |>
      ggplot(aes(.data[[pc_cols[1]]], .data[[pc_cols[2]]], fill = tissue)) +
      geom_point(shape = 21, size = 2.8, color = "black", stroke = 0.25, alpha = 0.9) +
      labs(
        title = "CIBERSORTx immune-composition principal components",
        x = pc_cols[1],
        y = pc_cols[2],
        fill = "Tissue"
      ) +
      theme_pub(10)
  }
  
  if (length(frac_cols) >= 2) {
    vars <- sort(vapply(dat[frac_cols], var, numeric(1), na.rm = TRUE), decreasing = TRUE)
    sel <- names(vars)[seq_len(min(10, length(vars)))]
    
    frac_long <- dat |>
      select(sample_id, tissue, all_of(sel)) |>
      pivot_longer(all_of(sel), names_to = "cell_fraction", values_to = "fraction") |>
      filter(is.finite(fraction)) |>
      mutate(cell_fraction = str_replace_all(cell_fraction, "\\.", " "))
    
    plots[[length(plots) + 1]] <- ggplot(frac_long, aes(tissue, fraction, fill = tissue)) +
      geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.82, linewidth = 0.35) +
      geom_jitter(width = 0.14, size = 0.55, alpha = 0.45) +
      facet_wrap(~ cell_fraction, scales = "free_y", ncol = 5) +
      labs(
        title = "Most variable CIBERSORTx-inferred immune fractions",
        x = NULL,
        y = "Relative fraction",
        fill = "Tissue"
      ) +
      theme_pub(8.5) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 8, face = "bold")
      )
    
    write_tsv(frac_long, "results/tables/suppfig_s5_cibersortx_fraction_plot_data.tsv")
  }
  
  if (length(plots) == 0) stop("No CIBERSORTx PC or fraction columns found.")
  
  p <- if (length(plots) == 1) plots[[1]] else combine_plots(plots[[1]], plots[[2]], ncol = 1)
  
  save_pub(
    p,
    "Supplementary_Figure_S5_CIBERSORTx_immune_composition",
    width = 10.5,
    height = ifelse(length(plots) == 1, 5.8, 10.5)
  )
  
  write_tsv(dat, "results/tables/suppfig_s5_cibersortx_input_data.tsv")
}

# ---------------------------
# Supplementary Figure S6
# TCGA Cox forest plot
# ---------------------------

make_s6 <- function() {
  cox_file <- "results/tables/tcga_lihc_survival_cox_models.tsv"
  if (!file.exists(cox_file)) stop("Missing: ", cox_file)
  
  cox <- read_tsv(cox_file, show_col_types = FALSE)
  
  term_col <- find_col(cox, c("^term$"))
  score_col <- find_col(cox, c("^score$"), required = FALSE)
  model_col <- find_col(cox, c("^model$"))
  hr_col <- find_col(cox, c("^estimate$", "^HR$", "hazard"))
  low_col <- find_col(cox, c("^conf.low$", "ci_low", "lower"))
  high_col <- find_col(cox, c("^conf.high$", "ci_high", "upper"))
  p_col <- find_col(cox, c("^p.value$", "^p_value$", "^p$"))
  
  plot_dat <- cox |>
    mutate(
      term = as.character(.data[[term_col]]),
      score = if (!is.na(score_col)) as.character(.data[[score_col]]) else term,
      model = as.character(.data[[model_col]]),
      HR = as.numeric(.data[[hr_col]]),
      conf_low = as.numeric(.data[[low_col]]),
      conf_high = as.numeric(.data[[high_col]]),
      p_value = as.numeric(.data[[p_col]])
    ) |>
    filter(
      str_detect(term, regex("ProlifHubScore|HepLossScore|HCCStateScore", TRUE)) |
        score %in% c("ProlifHubScore", "HepLossScore", "HCCStateScore")
    ) |>
    mutate(
      score = factor(score, levels = c("ProlifHubScore", "HepLossScore", "HCCStateScore")),
      model = recode(
        model,
        "score_only" = "Score only",
        "age_sex_adjusted" = "Age/sex adjusted",
        "age_sex_stage_adjusted" = "Age/sex/stage adjusted"
      ),
      model = factor(model, levels = c("Score only", "Age/sex adjusted", "Age/sex/stage adjusted"))
    ) |>
    filter(is.finite(HR), is.finite(conf_low), is.finite(conf_high))
  
  p <- ggplot(plot_dat, aes(HR, score, xmin = conf_low, xmax = conf_high)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey40") +
    geom_errorbarh(height = 0.18, linewidth = 0.45) +
    geom_point(size = 2.6) +
    scale_x_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
      labels = c("0.5", "0.75", "1", "1.5", "2", "3")
    ) +
    facet_wrap(~ model, ncol = 1) +
    labs(
      title = "TCGA-LIHC Cox model hazard ratios",
      x = "Hazard ratio per 1 SD increase in score",
      y = NULL
    ) +
    theme_pub(10) +
    theme(axis.text.y = element_text(face = "bold"))
  
  save_pub(p, "Supplementary_Figure_S6_TCGA_Cox_forest", width = 7.8, height = 7.0)
  
  write_tsv(plot_dat, "results/tables/suppfig_s6_tcga_cox_forest_plot_data.tsv")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "all"

run_one <- function(x) {
  message("\n--- Generating ", x, " ---")
  tryCatch(
    switch(
      x,
      S1 = make_s1(),
      S2 = make_s2(),
      S3 = make_s3(),
      S4 = make_s4(),
      S5 = make_s5(),
      S6 = make_s6(),
      stop("Unknown figure code: ", x)
    ),
    error = function(e) message("FAILED ", x, ": ", conditionMessage(e))
  )
}

if ("all" %in% args) {
  walk(c("S1", "S2", "S3", "S4", "S5", "S6"), run_one)
} else {
  walk(args, run_one)
}

message("\nDone.")