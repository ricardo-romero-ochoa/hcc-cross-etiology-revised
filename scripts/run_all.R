source("R/_shared.R")
ensure_dirs(); cfg <- read_config()


# Core GEO/discovery/validation pipeline.
# CIBERSORTx-adjusted regression is conditional because CIBERSORTx must be run externally.
steps <- c(
  "scripts/01_dataset_inventory.R",
  "scripts/02_download_and_curate_geo.R",
  "scripts/02c_audit_tissue_labels.R",
  "scripts/02b_audit_gene_symbol_overlap.R",

steps <- c(
  "scripts/01_dataset_inventory.R",
  "scripts/02_download_and_curate_geo.R",

  "scripts/03_discovery_hallmark_gsva.R",
  "scripts/03b_derive_hepatitis_axes.R",
  "scripts/04_meta_modules.R",
  "scripts/05_validate_geo_cohorts.R",
  "scripts/10_module_size_sensitivity.R",

  "scripts/08_cibersortx_export_GSE121248.R"
)

if (isTRUE(cfg$run_estimate)) {
  steps <- append(steps, "scripts/07_estimate_adjustment.R", after = length(steps))
}

if (isTRUE(cfg$run_tcga)) {
  steps <- append(steps, "scripts/06_tcga_lihc_pipeline.R", after = length(steps))
}

if (isTRUE(cfg$run_cibersortx_import)) {
  steps <- append(steps, "scripts/09_cibersortx_adjusted_regression_GSE121248.R", after = length(steps))
  steps <- append(steps, "scripts/11b_make_missing_figures.R", after = length(steps))
}

# Manuscript/revision outputs that can be regenerated from available intermediate files.
steps <- c(
  steps,
  "scripts/11_make_revision_figures.R",
  "scripts/12_make_manuscript_tables.R",
  "scripts/13_summarize_revision_results.R",
  "scripts/14_make_supplementary_figures.R",
  "scripts/15_make_supplementary_tables.R"
)

for (s in steps) {
  log_msg("Running", s)
  tryCatch(
    sys.source(s, envir = new.env(parent = globalenv())),
    error = function(e) {
      dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
      writeLines(conditionMessage(e), paste0("results/logs/ERROR_", basename(s), ".log"))
      stop(e)
    }
  )
}

message("Pipeline completed.")
message("If CIBERSORTx was not imported, run CIBERSORTx externally, set run_cibersortx_import: true, rerun the pipeline or run scripts/09_cibersortx_adjusted_regression_GSE121248.R manually.")
message("Run the final repository audit with: Rscript scripts/16_revision_repository_audit.R --strict")

  "scripts/08_cibersortx_export_GSE121248.R",
  "scripts/11_make_revision_figures.R",
  "scripts/12_make_manuscript_tables.R"
)
if (isTRUE(cfg$run_estimate)) steps <- append(steps, "scripts/07_estimate_adjustment.R", after = 8)
if (isTRUE(cfg$run_tcga)) steps <- append(steps, "scripts/06_tcga_lihc_pipeline.R", after = 8)
for (s in steps) {
  log_msg("Running", s)
  tryCatch(sys.source(s, envir = new.env(parent = globalenv())), error = function(e) {
    writeLines(conditionMessage(e), paste0("results/logs/ERROR_", basename(s), ".log"))
    stop(e)
  })
}
message("Pipeline completed. Run scripts/09_cibersortx_adjusted_regression_GSE121248.R after generating the external CIBERSORTx output file.")

