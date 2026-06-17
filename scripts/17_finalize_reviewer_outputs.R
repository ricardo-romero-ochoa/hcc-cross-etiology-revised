#!/usr/bin/env Rscript
# =============================================================================
# 17_finalize_reviewer_outputs.R
#
# Purpose:
#   Hot-fix/finalization script for the major-revision repository. Use this after
#   the core pipeline has completed to generate reviewer-facing audit tables,
#   manuscript-table aliases, supplementary-table index, and then run the strict
#   repository audit.
#
# Usage:
#   Rscript scripts/17_finalize_reviewer_outputs.R
# =============================================================================

run_script <- function(path) {
  message("[finalize] Running ", path)
  if (!file.exists(path)) stop("Missing script: ", path, call. = FALSE)
  sys.source(path, envir = new.env(parent = globalenv()))
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/supplementary_tables", recursive = TRUE, showWarnings = FALSE)

# Audits not originally included in run_all.R in earlier corrected builds.
run_script("scripts/02c_audit_tissue_labels.R")
run_script("scripts/02b_audit_gene_symbol_overlap.R")

# Regenerate manuscript and supplementary table outputs/aliases.
run_script("scripts/12_make_manuscript_tables.R")
run_script("scripts/15_make_supplementary_tables.R")

# Defensive alias creation for older script versions.
aliases <- list(
  c("manuscript/tables/Table1_meta_analysis_summary.csv", "manuscript/tables/Table1_meta_summary.csv"),
  c("manuscript/tables/Table2_HBV_injury_adjusted_models.csv", "manuscript/tables/Table2_HBV_injury_regression.csv"),
  c("manuscript/tables/Table3_TCGA_LIHC_Cox_models.csv", "manuscript/tables/Table3_TCGA_Cox.csv")
)
for (pair in aliases) {
  src <- pair[[1]]; dst <- pair[[2]]
  if (file.exists(src) && !file.exists(dst)) {
    file.copy(src, dst, overwrite = TRUE)
    message("[finalize] Created alias: ", dst)
  }
}

if (!file.exists("manuscript/supplementary_tables/Supplementary_Table_Index.tsv") &&
    file.exists("results/tables/supplementary_tables_file_status.tsv")) {
  status <- readr::read_tsv("results/tables/supplementary_tables_file_status.tsv", show_col_types = FALSE)
  readr::write_tsv(status, "manuscript/supplementary_tables/Supplementary_Table_Index.tsv")
  readr::write_csv(status, "manuscript/supplementary_tables/Supplementary_Table_Index.csv")
  message("[finalize] Created supplementary table index from status table.")
}

# Strict audit.
message("[finalize] Running strict repository audit")
source("scripts/16_revision_repository_audit.R", local = new.env(parent = globalenv()))
message("[finalize] Done. If the strict audit still fails, inspect results/revision_audit/repository_manifest_audit.tsv")
