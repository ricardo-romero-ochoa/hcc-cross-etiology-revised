# =============================================================================
# Repository audit for INSILI-D-26-00441 major revision
#
# Purpose:
#   Verifies that reviewer-requested reproducibility objects are discoverable in
#   the public repository/release. The audit reads REPOSITORY_MANIFEST.csv and
#   writes machine-readable audit outputs to results/revision_audit/.
#
# Usage:
#   Source-only audit before running generated outputs:
#     Rscript scripts/16_revision_repository_audit.R
#
#   Final strict audit before GitHub release / Zenodo version:
#     Rscript scripts/16_revision_repository_audit.R --strict
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
strict <- any(args %in% c("--strict", "strict", "TRUE", "true"))
manifest_path <- "REPOSITORY_MANIFEST.csv"
out_dir <- "results/revision_audit"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

if (!file.exists(manifest_path)) {
  fail("Missing ", manifest_path, ". Add the reviewer-facing repository manifest before running this audit.")
}

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("category", "file_or_folder", "description", "required_for_reproducibility", "object_type", "strict_audit")
missing_cols <- setdiff(required_cols, names(manifest))
if (length(missing_cols) > 0) {
  fail("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

is_yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1", "required")

manifest$required_flag <- is_yes(manifest$required_for_reproducibility)
manifest$strict_flag <- is_yes(manifest$strict_audit)
manifest$exists <- file.exists(manifest$file_or_folder)
manifest$size_bytes <- vapply(manifest$file_or_folder, function(p) {
  if (!file.exists(p)) return(NA_real_)
  if (dir.exists(p)) return(NA_real_)
  as.numeric(file.info(p)$size)
}, numeric(1))
manifest$nonempty <- ifelse(is.na(manifest$size_bytes), manifest$exists, manifest$size_bytes > 0)
manifest$audit_mode <- if (strict) "strict" else "source_only"
manifest$checked_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# In source-only mode, require source/config/manifest/audit objects only.
# In strict mode, require all entries marked strict_audit == yes.
if (strict) {
  audit_df <- manifest[manifest$required_flag & manifest$strict_flag, , drop = FALSE]
} else {
  audit_df <- manifest[manifest$required_flag & manifest$object_type %in% c("source", "config", "audit"), , drop = FALSE]
}

missing <- audit_df[!audit_df$nonempty, , drop = FALSE]

utils::write.csv(manifest, file.path(out_dir, "repository_manifest_audit.csv"), row.names = FALSE)
utils::write.table(manifest, file.path(out_dir, "repository_manifest_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

summary_df <- data.frame(
  audit_mode = if (strict) "strict" else "source_only",
  checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  manifest_rows = nrow(manifest),
  checked_rows = nrow(audit_df),
  missing_rows = nrow(missing),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_df, file.path(out_dir, "repository_manifest_audit_summary.csv"), row.names = FALSE)

if (nrow(missing) > 0) {
  message("Repository audit failed. Missing required objects:")
  print(missing[, c("category", "file_or_folder", "description", "object_type")], row.names = FALSE)
  if (strict) {
    message("\nStrict audit should be run after the full pipeline has generated reviewer-facing outputs.")
    message("The TCGA Cox output expected by the reviewer/manuscript is: results/tables/tcga_lihc_survival_cox_models.tsv")
  }
  fail("Repository audit failed with ", nrow(missing), " missing required object(s).")
}

message("Repository audit passed in ", if (strict) "strict" else "source-only", " mode.")
message("Audit outputs written to ", out_dir, "/")
