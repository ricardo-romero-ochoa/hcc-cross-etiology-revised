#!/usr/bin/env Rscript

# Generate revised supplementary tables S2-S11 for the HCC cross-etiology manuscript.
# S1 is assumed to be the merged dataset inventory/metadata table prepared separately.
#
# Outputs:
#   manuscript/supplementary_tables/Supplementary_Tables_S2_to_S11.docx
#   manuscript/supplementary_tables/Supplementary_Table_S*.docx
#   manuscript/supplementary_tables/Supplementary_Table_S*.tsv
#   results/tables/supplementary_tables_file_status.tsv
#
# Required packages:
#   install.packages(c("officer", "flextable", "dplyr", "readr", "tibble", "stringr", "purrr", "tidyr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(officer)
  library(flextable)
})

if (file.exists("R/_shared.R")) source("R/_shared.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

out_dir <- "manuscript/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# General helpers
# -----------------------------

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
    if (required) stop("Could not find ", label, ". Available columns: ", paste(colnames(df), collapse = ", "))
    return(NA_character_)
  }
  hits[1]
}

safe_read_tsv <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

clean_names_for_word <- function(x) {
  x <- as.data.frame(x)
  colnames(x) <- colnames(x) |>
    str_replace_all("_", " ") |>
    str_replace_all("\\.", " ") |>
    str_squish()
  x
}

format_for_word <- function(x) {
  x <- x |> as.data.frame()
  x[] <- lapply(x, function(col) {
    if (is.numeric(col)) {
      out <- ifelse(
        is.na(col),
        NA_character_,
        ifelse(abs(col) > 0 & (abs(col) < 0.001 | abs(col) >= 10000),
               formatC(col, format = "e", digits = 3),
               formatC(col, format = "f", digits = 4))
      )
      return(out)
    }
    as.character(col)
  })
  clean_names_for_word(x)
}

write_tsv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df, path)
  invisible(path)
}

add_ft_to_doc <- function(doc, df, title, caption, source_note = NULL, max_rows = 80) {
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, caption, style = "Normal")
  if (!is.null(source_note)) doc <- body_add_par(doc, source_note, style = "Normal")

  n_total <- nrow(df)
  shown <- df
  if (n_total > max_rows) shown <- df |> slice_head(n = max_rows)

  if (nrow(shown) == 0) {
    doc <- body_add_par(doc, "No rows available.", style = "Normal")
  } else {
    ft <- flextable(format_for_word(shown)) |>
      theme_booktabs() |>
      fontsize(size = 7, part = "all") |>
      fontsize(size = 7, part = "header") |>
      align(align = "left", part = "all") |>
      valign(valign = "top", part = "all") |>
      autofit()
    ft <- set_table_properties(ft, layout = "autofit", width = 1)
    doc <- body_add_flextable(doc, ft)
  }
  if (n_total > max_rows) {
    doc <- body_add_par(doc, paste0("Preview shown: first ", max_rows, " of ", n_total, " rows. The full table is provided as a TSV file in the same folder."), style = "Normal")
  }
  doc <- body_add_par(doc, "", style = "Normal")
  doc
}

write_individual_docx <- function(df, table_id, title, caption, source_note = NULL, max_rows = 80) {
  doc <- read_docx()
  doc <- add_ft_to_doc(doc, df, paste0("Supplementary Table ", table_id, ". ", title), caption, source_note, max_rows)
  path <- file.path(out_dir, paste0("Supplementary_Table_", table_id, ".docx"))
  print(doc, target = path)
  invisible(path)
}

add_missing_table <- function(table_id, title, missing_sources, instructions) {
  tibble(
    table_id = table_id,
    title = title,
    status = "missing_source_files",
    missing_sources = paste(missing_sources, collapse = "; "),
    instructions = instructions
  )
}

status_rows <- list()
tables <- list()

register_table <- function(table_id, title, caption, df, full_tsv_path, source_files, max_rows_word = 80) {
  write_tsv_safe(df, full_tsv_path)
  write_individual_docx(df, table_id, title, caption,
                        source_note = paste("Source file(s):", paste(source_files, collapse = "; ")),
                        max_rows = max_rows_word)
  tables[[length(tables) + 1]] <<- list(
    table_id = table_id,
    title = title,
    caption = caption,
    df = df,
    full_tsv_path = full_tsv_path,
    source_files = source_files,
    max_rows_word = max_rows_word
  )
  status_rows[[length(status_rows) + 1]] <<- tibble(
    table_id = table_id,
    title = title,
    status = "generated",
    rows = nrow(df),
    columns = ncol(df),
    full_tsv_path = full_tsv_path,
    docx_path = file.path(out_dir, paste0("Supplementary_Table_", table_id, ".docx")),
    source_files = paste(source_files, collapse = "; "),
    instructions = NA_character_
  )
}

# -----------------------------
# S2. Complete Hallmark GSVA differential activity
# -----------------------------

make_s2 <- function() {
  combined <- first_existing(c(
    "results/tables/discovery_hallmark_combined.tsv",
    "results/tables/hallmark_discovery_combined.tsv"
  ))

  if (!is.na(combined)) {
    df <- read_tsv(combined, show_col_types = FALSE)
  } else {
    files <- c(
      "results/tables/GSE121248_hallmark_limma.tsv",
      "results/tables/GSE41804_hallmark_limma.tsv"
    )
    if (!all(file.exists(files))) {
      status_rows[[length(status_rows) + 1]] <<- add_missing_table(
        "S2", "Complete Hallmark GSVA differential activity",
        c(combined, files[!file.exists(files)]),
        "Run: Rscript scripts/03_discovery_hallmark_gsva.R"
      )
      return(invisible(NULL))
    }
    df <- map_dfr(files, function(f) {
      acc <- str_extract(basename(f), "GSE[0-9]+")
      read_tsv(f, show_col_types = FALSE) |> mutate(dataset = acc, .before = 1)
    })
  }

  set_col <- find_col(df, c("Hallmark", "hallmark", "gs_name", "pathway", "gene_set"), required = FALSE)
  logfc_col <- find_col(df, c("^logFC$", "log_fc", "estimate"), required = FALSE)
  fdr_col <- find_col(df, c("^FDR$", "adj.P.Val", "adj_p", "padj"), required = FALSE)
  p_col <- find_col(df, c("^P.Value$", "^p_value$", "^p$", "pval"), required = FALSE)

  out <- df
  if (!is.na(set_col)) out <- out |> rename(hallmark_gene_set = all_of(set_col))
  if (!is.na(logfc_col) && logfc_col != "logFC") out <- out |> rename(logFC = all_of(logfc_col))
  if (!is.na(p_col) && p_col != "p_value") out <- out |> rename(p_value = all_of(p_col))
  if (!is.na(fdr_col) && fdr_col != "FDR") out <- out |> rename(FDR = all_of(fdr_col))

  register_table(
    "S2",
    "Complete Hallmark GSVA differential activity in discovery cohorts",
    "Full Hallmark pathway-level tumor versus adjacent/non-tumor differential activity results for the HBV-HCC and HCV-HCC discovery cohorts.",
    out,
    file.path(out_dir, "Supplementary_Table_S2_Hallmark_GSVA_results.tsv"),
    source_files = if (!is.na(combined)) combined else files,
    max_rows_word = 120
  )
}

# -----------------------------
# S3. Full meta-analysis with heterogeneity
# -----------------------------

make_s3 <- function() {
  f <- "results/tables/meta_analysis_with_heterogeneity.tsv"
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S3", "Full gene-level meta-analysis with heterogeneity",
      f,
      "Run: Rscript scripts/04_meta_modules.R"
    )
    return(invisible(NULL))
  }

  preferred_cols <- c(
    "gene", "meta_status",
    "logFC_GSE121248", "SE_GSE121248", "FDR_GSE121248",
    "logFC_GSE41804", "SE_GSE41804", "FDR_GSE41804",
    "meta_logFC_FE", "meta_se_FE", "meta_z_FE", "meta_p_FE", "meta_FDR_FE",
    "Q", "Q_p", "I2", "tau2_REML", "meta_logFC_RE", "meta_se_RE", "meta_p_RE", "meta_FDR_RE",
    "direction_concordant"
  )
  out <- df |> select(any_of(preferred_cols), everything())

  register_table(
    "S3",
    "Full gene-level meta-analysis with heterogeneity statistics",
    "Complete gene-level fixed-effect meta-analysis and random-effects heterogeneity output across GSE121248 and GSE41804. The Word file contains a preview because the full table is large; the complete table is provided as TSV.",
    out,
    file.path(out_dir, "Supplementary_Table_S3_meta_analysis_with_heterogeneity.tsv"),
    source_files = f,
    max_rows_word = 80
  )
}

# -----------------------------
# S4. ProlifHub and HepLoss module definitions
# -----------------------------

make_s4 <- function() {
  f <- first_existing(c(
    "results/tables/module_gene_panels_top20.tsv",
    "results/tables/module_gene_panels_top20.tsv"
  ))
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S4", "ProlifHub and HepLoss module definitions",
      f,
      "Run: Rscript scripts/04_meta_modules.R"
    )
    return(invisible(NULL))
  }

  meta_file <- "results/tables/meta_analysis_with_heterogeneity.tsv"
  meta <- safe_read_tsv(meta_file)

  gene_col <- find_col(df, c("^gene$", "symbol", "gene_symbol"), required = FALSE)
  if (!is.na(gene_col) && !is.null(meta) && "gene" %in% colnames(meta)) {
    out <- df |>
      rename(gene = all_of(gene_col)) |>
      left_join(meta |> select(any_of(c("gene", "meta_logFC_FE", "meta_FDR_FE", "I2", "tau2_REML"))), by = "gene")
  } else {
    out <- df
  }

  register_table(
    "S4",
    "ProlifHub and HepLoss module definitions",
    "Genes included in the compact ProlifHub and HepLoss modules, with available meta-analysis statistics where present.",
    out,
    file.path(out_dir, "Supplementary_Table_S4_module_definitions.tsv"),
    source_files = c(f, if (!is.null(meta)) meta_file else NA_character_) |> na.omit(),
    max_rows_word = 120
  )
}

# -----------------------------
# S5. Multi-cohort validation statistics
# -----------------------------

make_s5 <- function() {
  f <- first_existing(c(
    "results/tables/geo_validation_summary.tsv",
    "manuscript/tables/SuppTable_multicohort_validation.csv"
  ))
  if (is.na(f)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S5", "Multi-cohort module-score validation statistics",
      "results/tables/geo_validation_summary.tsv",
      "Run: Rscript scripts/05_validate_geo_cohorts.R"
    )
    return(invisible(NULL))
  }

  df <- if (str_ends(f, ".csv")) read_csv(f, show_col_types = FALSE) else read_tsv(f, show_col_types = FALSE)

  register_table(
    "S5",
    "Multi-cohort module-score validation statistics",
    "Tumor versus non-tumor validation statistics for ProlifHubScore, HepLossScore, and HCCStateScore across independent GEO HCC cohorts.",
    df,
    file.path(out_dir, "Supplementary_Table_S5_multicohort_validation.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# S6. Module-size robustness statistics
# -----------------------------

make_s6 <- function() {
  f <- first_existing(c(
    "results/tables/module_size_robustness.tsv",
    "results/tables/module_size_robustness_summary.tsv",
    "results/tables/module_size_sensitivity.tsv"
  ))
  if (is.na(f)) {
    hits <- list.files("results/tables", pattern = "module.*size|size.*robust|module.*sensitivity", full.names = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) f <- hits[1]
  }
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S6", "Module-size robustness statistics",
      "results/tables/module_size_robustness.tsv",
      "Run: Rscript scripts/10_module_size_sensitivity.R"
    )
    return(invisible(NULL))
  }

  register_table(
    "S6",
    "Module-size robustness statistics",
    "HCCStateScore validation results across alternative module sizes used to assess robustness to the top-20 module definition.",
    df,
    file.path(out_dir, "Supplementary_Table_S6_module_size_robustness.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# S7. HBV_INJURY derivation table
# -----------------------------

make_s7 <- function() {
  f <- "results/tables/GSE83148_HBV_INJURY_derivation_full.tsv"
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S7", "HBV_INJURY derivation from GSE83148",
      f,
      "Run: Rscript scripts/03b_derive_hepatitis_axes.R, or the equivalent HBV injury derivation script."
    )
    return(invisible(NULL))
  }

  # Add top-N membership flags if absent.
  if (all(c("gene", "beta_injury_index", "FDR") %in% colnames(df))) {
    stat_col <- find_col(df, c("^t$", "stat"), required = FALSE)
    if (is.na(stat_col)) stat_col <- "beta_injury_index"
    ranked <- df |>
      filter(beta_injury_index > 0, FDR < 0.10) |>
      mutate(.ranking_stat = abs(.data[[stat_col]])) |>
      arrange(FDR, desc(.ranking_stat), desc(beta_injury_index)) |>
      mutate(hbv_injury_rank = row_number()) |>
      select(gene, hbv_injury_rank)
    df <- df |>
      left_join(ranked, by = "gene") |>
      mutate(
        included_TOP_200 = !is.na(hbv_injury_rank) & hbv_injury_rank <= 200,
        included_TOP_500 = !is.na(hbv_injury_rank) & hbv_injury_rank <= 500,
        included_TOP_1000 = !is.na(hbv_injury_rank) & hbv_injury_rank <= 1000,
        included_TOP_2000 = !is.na(hbv_injury_rank) & hbv_injury_rank <= 2000,
        included_TOP_5000 = !is.na(hbv_injury_rank) & hbv_injury_rank <= 5000,
        included_EXTENDED_7792 = !is.na(hbv_injury_rank)
      ) |>
      arrange(hbv_injury_rank, FDR)
  }

  register_table(
    "S7",
    "HBV_INJURY derivation from GSE83148",
    "Gene-level association with the continuous ALT/AST/HBV-DNA-derived HBV_INJURY_INDEX in GSE83148, including ranked top-N membership indicators where available. The Word file contains a preview; the complete table is provided as TSV.",
    df,
    file.path(out_dir, "Supplementary_Table_S7_HBV_INJURY_derivation.tsv"),
    source_files = f,
    max_rows_word = 80
  )
}

# -----------------------------
# S8. HBV_INJURY top-N/CIBERSORTx regression
# -----------------------------

make_s8 <- function() {
  f <- first_existing(c(
    "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv",
    "results/tables/gse121248_hbv_injury_topN_extended_regression_tissue.tsv"
  ))
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S8", "HBV_INJURY top-N and CIBERSORTx-adjusted regression",
      "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv",
      "Run CIBERSORTx externally, then run: Rscript scripts/09_cibersortx_adjusted_regression_GSE121248.R"
    )
    return(invisible(NULL))
  }

  register_table(
    "S8",
    "HBV_INJURY top-N and CIBERSORTx-adjusted regression",
    "Nested regression estimates for HBV_INJURY gene-set definitions in GSE121248 before and after proliferation and CIBERSORTx-inferred immune-composition adjustment.",
    df,
    file.path(out_dir, "Supplementary_Table_S8_HBV_INJURY_topN_CIBERSORTx_regression.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# S9. CIBERSORTx fraction summaries
# -----------------------------

make_s9 <- function() {
  f <- first_existing(c(
    "results/tables/gse121248_cibersortx_cell_fraction_summary.tsv",
    "results/tables/suppfig_s5_cibersortx_fraction_plot_data.tsv"
  ))
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S9", "CIBERSORTx inferred immune-fraction summaries",
      "results/tables/gse121248_cibersortx_cell_fraction_summary.tsv",
      "Run CIBERSORTx externally, then run: Rscript scripts/09_cibersortx_adjusted_regression_GSE121248.R"
    )
    return(invisible(NULL))
  }

  register_table(
    "S9",
    "CIBERSORTx inferred immune-fraction summaries",
    "Summary statistics for CIBERSORTx-inferred immune-cell fractions by tissue class in GSE121248.",
    df,
    file.path(out_dir, "Supplementary_Table_S9_CIBERSORTx_fraction_summaries.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# S10. TCGA Cox model full output
# -----------------------------

make_s10 <- function() {
  f <- "results/tables/tcga_lihc_survival_cox_models.tsv"
  df <- safe_read_tsv(f)
  if (is.null(df)) {
    status_rows[[length(status_rows) + 1]] <<- add_missing_table(
      "S10", "TCGA-LIHC Cox model full output",
      f,
      "Run: Rscript scripts/06_tcga_lihc_pipeline.R, or the patched TCGA-LIHC Cox script."
    )
    return(invisible(NULL))
  }

  register_table(
    "S10",
    "TCGA-LIHC Cox model full output",
    "Full Cox proportional-hazards model output for TCGA-LIHC module scores, including score-only, age/sex-adjusted, and age/sex/stage-adjusted models.",
    df,
    file.path(out_dir, "Supplementary_Table_S10_TCGA_LIHC_Cox_models.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# S11. Session information
# -----------------------------

make_s11 <- function() {
  f <- first_existing(c(
    "results/tables/session_info.tsv",
    "results/tables/software_session_info.tsv"
  ))

  if (!is.na(f)) {
    df <- read_tsv(f, show_col_types = FALSE)
  } else {
    si <- utils::sessionInfo()
    pkg_df <- as.data.frame(si$otherPkgs %||% list())
    # Build simpler package-version table from installed packages used by pipeline.
    pkgs <- c("limma", "GSVA", "msigdbr", "metafor", "pROC", "TCGAbiolinks", "survival", "broom", "tidyverse", "dplyr", "readr", "ggplot2", "officer", "flextable")
    ip <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
    df <- tibble(
      item = c("R.version", "platform", "running", "date", paste0("package:", pkgs)),
      value = c(
        R.version.string,
        si$platform,
        si$running,
        as.character(Sys.Date()),
        ip$Version[match(pkgs, ip$Package)]
      )
    )
    write_tsv_safe(df, "results/tables/session_info.tsv")
    f <- "results/tables/session_info.tsv"
  }

  register_table(
    "S11",
    "Software and session information",
    "R, platform, date, and package-version information for reproducibility.",
    df,
    file.path(out_dir, "Supplementary_Table_S11_session_information.tsv"),
    source_files = f,
    max_rows_word = 120
  )
}

# -----------------------------
# Run all table builders
# -----------------------------

builders <- list(make_s2, make_s3, make_s4, make_s5, make_s6, make_s7, make_s8, make_s9, make_s10, make_s11)
walk(builders, function(fun) {
  tryCatch(fun(), error = function(e) {
    warning("Table generation failed: ", conditionMessage(e))
  })
})

# -----------------------------
# Combined Word document
# -----------------------------

if (length(tables) > 0) {
  doc <- read_docx()
  doc <- body_add_par(doc, "Supplementary Tables S2-S11", style = "heading 1")
  doc <- body_add_par(doc, "S1 is assumed to be the merged dataset inventory and curated metadata table prepared separately. Large supplementary tables are represented in this Word document as previews; complete TSV versions are written in manuscript/supplementary_tables/.", style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")

  for (tbl in tables) {
    doc <- add_ft_to_doc(
      doc,
      tbl$df,
      paste0("Supplementary Table ", tbl$table_id, ". ", tbl$title),
      tbl$caption,
      source_note = paste0("Complete TSV: ", tbl$full_tsv_path, ". Source file(s): ", paste(tbl$source_files, collapse = "; ")),
      max_rows = tbl$max_rows_word
    )
  }

  combined_docx <- file.path(out_dir, "Supplementary_Tables_S2_to_S11.docx")
  print(doc, target = combined_docx)
  message("Wrote combined DOCX: ", combined_docx)
}

status <- bind_rows(status_rows)
write_tsv_safe(status, "results/tables/supplementary_tables_file_status.tsv")
write_tsv_safe(status, file.path(out_dir, "supplementary_tables_file_status.tsv"))

# Reviewer-facing supplementary-table index required by the repository manifest.
# This is intentionally separate from the status table so reviewers can quickly
# locate each supplementary table and identify missing external-source outputs.
supp_index <- status |>
  dplyr::transmute(
    table_id = table_id,
    title = title,
    status = status,
    rows = dplyr::if_else(is.na(rows), NA_integer_, as.integer(rows)),
    columns = dplyr::if_else(is.na(columns), NA_integer_, as.integer(columns)),
    tsv_path = full_tsv_path,
    docx_path = docx_path,
    source_files = source_files,
    instructions = instructions
  )
write_tsv_safe(supp_index, file.path(out_dir, "Supplementary_Table_Index.tsv"))
readr::write_csv(supp_index, file.path(out_dir, "Supplementary_Table_Index.csv"))

message("Supplementary table status:")
print(status |> select(any_of(c("table_id", "title", "status", "rows", "instructions"))))

message("Done.")
