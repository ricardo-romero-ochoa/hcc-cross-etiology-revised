# =============================================================================
# Minor revision reviewer-facing outputs
#
# This script should be run after the full pipeline has completed. It creates
# reviewer-facing tables requested in the minor revision and adds
# proportional-hazards diagnostics for TCGA-LIHC Cox models using cox.zph().
# =============================================================================

source("R/_shared.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(survival)
  library(purrr)
})
ensure_dirs()
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/supplementary_tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Dataset summary table with required reviewer fields
# -----------------------------
platform_map <- c(
  GSE121248="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  GSE41804="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  GSE83148="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  GSE38941="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  GSE14520="GPL3921 + GPL571 (Affymetrix HT HG-U133A / HG-U133A 2.0)",
  GSE25097="GPL10687 (Rosetta/Merck Human RSTA Custom Affymetrix 1.0)",
  GSE36376="GPL10558 (Illumina HumanHT-12 V4.0)",
  GSE76427="GPL10558 (Illumina HumanHT-12 V4.0)",
  GSE57957="GPL10558 (Illumina HumanHT-12 V4.0)",
  GSE45267="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  GSE112790="GPL570 (Affymetrix HG-U133 Plus 2.0)",
  `TCGA-LIHC`="RNA-seq STAR-counts / GDC-TCGA"
)

tissue_type_map <- c(
  GSE121248="HCC tumor and paired adjacent/non-tumor liver",
  GSE41804="HCC tumor and adjacent/non-tumor liver",
  GSE83148="Chronic HBV hepatitis liver and healthy liver; no HCC tumors",
  GSE38941="HBV-associated acute liver failure explanted liver and normal donor liver; no HCC tumors",
  GSE14520="HCC tumor and adjacent/non-tumor liver",
  GSE25097="HCC tumor and non-tumor/adjacent/healthy/cirrhotic liver where labels were recoverable",
  GSE36376="HCC tumor and non-tumor/adjacent liver",
  GSE76427="HCC tumor and adjacent/non-tumor liver",
  GSE57957="HCC tumor and adjacent/non-tumor liver",
  GSE45267="HCC tumor and adjacent/non-tumor liver",
  GSE112790="HCC tumor and limited adjacent/non-tumor liver",
  `TCGA-LIHC`="Primary HCC tumor RNA-seq with clinical follow-up; solid normal samples parsed when available"
)

inventory <- readr::read_tsv("results/tables/dataset_inventory.tsv", show_col_types = FALSE)
label_files <- list.files("results/tables", pattern = "^curated_metadata_.*\\.tsv$", full.names = TRUE)
label_counts <- purrr::map_dfr(label_files, function(f) {
  acc <- stringr::str_match(basename(f), "curated_metadata_(.*)\\.tsv")[,2]
  dat <- readr::read_tsv(f, show_col_types = FALSE)
  tissue_col <- intersect(c("tissue", "tissue_class", "label"), names(dat))[1]
  if (is.na(tissue_col)) return(tibble(accession = acc, n_tumor = NA_integer_, n_non_tumor = NA_integer_))
  tibble(
    accession = acc,
    n_tumor = sum(dat[[tissue_col]] == "tumor", na.rm = TRUE),
    n_non_tumor = sum(dat[[tissue_col]] %in% c("non_tumor", "adjacent", "normal"), na.rm = TRUE)
  )
})

cohort_summary <- inventory |>
  left_join(label_counts, by = "accession") |>
  mutate(
    platform = unname(platform_map[accession]),
    tissue_type = unname(tissue_type_map[accession]),
    inclusion_exclusion_rationale = case_when(
      accession == "GSE83148" ~ "Included for HBV injury-axis derivation only; excluded from tumor/non-tumor HCC validation.",
      accession == "GSE38941" ~ "Included only as contextual HBV liver-injury reference; excluded from HCC tumor/non-tumor validation.",
      TRUE ~ paste0("Included: ", inclusion_reason)
    ),
    use_in_analysis = role
  ) |>
  select(cohort_accession = accession, etiology, platform, tissue_type, n_tumor, n_non_tumor, use_in_analysis, inclusion_exclusion_rationale)

if (file.exists("results/tables/tcga_lihc_survival_model_input.tsv")) {
  tcga_surv <- readr::read_tsv("results/tables/tcga_lihc_survival_model_input.tsv", show_col_types = FALSE)
  cohort_summary <- bind_rows(cohort_summary, tibble(
    cohort_accession = "TCGA-LIHC",
    etiology = "mixed/clinical TCGA-LIHC",
    platform = platform_map[["TCGA-LIHC"]],
    tissue_type = tissue_type_map[["TCGA-LIHC"]],
    n_tumor = nrow(tcga_surv),
    n_non_tumor = NA_integer_,
    use_in_analysis = "RNA-seq module scoring and overall-survival Cox modeling",
    inclusion_exclusion_rationale = "Included for independent RNA-seq prognostic-coherence assessment and clinical modeling; not used to derive GEO microarray modules."
  ))
}
readr::write_csv(cohort_summary, "manuscript/supplementary_tables/Supplementary_Table_S1_dataset_summary_required_fields.csv")
readr::write_tsv(cohort_summary, "manuscript/supplementary_tables/Supplementary_Table_S1_dataset_summary_required_fields.tsv")

# -----------------------------
# Main validation table: HCCStateScore only
# -----------------------------
val <- readr::read_tsv("results/tables/geo_validation_summary.tsv", show_col_types = FALSE)
name_map <- names(val)
# tolerate either old or clean column names
get_col <- function(candidates) intersect(candidates, names(val))[1]
score_col <- get_col(c("score"))
dataset_col <- get_col(c("dataset", "accession"))
nt_col <- get_col(c("n_tumor", "n tumor"))
nn_col <- get_col(c("n_non_tumor", "n non tumor"))
delta_col <- get_col(c("delta_tumor_minus_nontumor", "delta tumor minus nontumor"))
p_col <- get_col(c("p_value", "p value", "p"))
fdr_col <- get_col(c("FDR", "fdr"))
auc_col <- get_col(c("AUC", "auc"))
auc_low_col <- get_col(c("AUC_low", "AUC low", "auc_low"))
auc_high_col <- get_col(c("AUC_high", "AUC high", "auc_high"))

main_val <- val |>
  filter(.data[[score_col]] == "HCCStateScore") |>
  transmute(
    cohort = .data[[dataset_col]],
    n_tumor = .data[[nt_col]],
    n_non_tumor = .data[[nn_col]],
    total_n = .data[[nt_col]] + .data[[nn_col]],
    tumor_non_tumor_delta = .data[[delta_col]],
    AUC = .data[[auc_col]],
    AUC_CI_low = .data[[auc_low_col]],
    AUC_CI_high = .data[[auc_high_col]],
    p_value = .data[[p_col]],
    FDR = .data[[fdr_col]]
  )
readr::write_csv(main_val, "manuscript/tables/Table2_validation_results.csv")
readr::write_tsv(main_val, "manuscript/supplementary_tables/Supplementary_Table_S5_main_validation_results.tsv")

# -----------------------------
# TCGA Cox models with proportional-hazards diagnostics
# -----------------------------
cox <- readr::read_tsv("results/tables/tcga_lihc_survival_cox_models.tsv", show_col_types = FALSE)
survdat_path <- "results/tables/tcga_lihc_survival_model_input.tsv"
if (!file.exists(survdat_path)) stop("Missing ", survdat_path, ". Re-run scripts/06_tcga_lihc_pipeline.R.")
survdat <- readr::read_tsv(survdat_path, show_col_types = FALSE) |>
  mutate(gender_tmp = as.factor(gender_tmp), stage_collapsed = factor(stage_collapsed, levels = c("I", "II", "III", "IV")))

fit_ph_one <- function(score_name, covars = character(), model_name = "score_only") {
  rhs <- paste0("scale(", score_name, ")")
  if (length(covars) > 0) rhs <- paste(rhs, paste(covars, collapse = " + "), sep = " + ")
  f <- as.formula(paste("survival::Surv(time_months, event) ~", rhs))
  dat <- survdat |> filter(!is.na(.data[[score_name]]), is.finite(time_months), time_months > 0)
  fit <- survival::coxph(f, data = dat)
  z <- survival::cox.zph(fit)
  zz <- as.data.frame(z$table) |> rownames_to_column("term")
  score_term <- zz |> filter(stringr::str_detect(term, paste0("scale\\(", score_name, "\\)")))
  global_term <- zz |> filter(term == "GLOBAL")
  tibble(
    score = score_name,
    model = model_name,
    ph_score_term_p = ifelse(nrow(score_term) == 1, score_term$p, NA_real_),
    ph_global_p = ifelse(nrow(global_term) == 1, global_term$p, NA_real_)
  )
}

ph <- bind_rows(
  lapply(c("ProlifHubScore", "HepLossScore", "HCCStateScore"), fit_ph_one, covars = character(), model_name = "score_only"),
  lapply(c("ProlifHubScore", "HepLossScore", "HCCStateScore"), fit_ph_one, covars = c("age_tmp", "gender_tmp"), model_name = "age_sex_adjusted"),
  lapply(c("ProlifHubScore", "HepLossScore", "HCCStateScore"), fit_ph_one, covars = c("age_tmp", "gender_tmp", "stage_collapsed"), model_name = "age_sex_stage_adjusted")
)
readr::write_tsv(ph, "results/tables/tcga_lihc_survival_ph_assumption.tsv")

covariate_map <- c(
  score_only = "None; score only",
  age_sex_adjusted = "Age and sex",
  age_sex_stage_adjusted = "Age, sex, and pathologic stage"
)
cox_main <- cox |>
  filter(stringr::str_detect(term, "^scale\\(")) |>
  left_join(ph, by = c("score", "model")) |>
  mutate(
    covariates = unname(covariate_map[model]),
    HR_per_SD = estimate,
    CI_95 = paste0(signif(conf.low, 3), "-", signif(conf.high, 3))
  ) |>
  select(score, model, covariates, n, events, HR_per_SD, conf.low, conf.high, p.value, ph_score_term_p, ph_global_p, formula)
readr::write_csv(cox_main, "manuscript/tables/Table4_TCGA_Cox_complete.csv")
readr::write_tsv(cox_main, "manuscript/supplementary_tables/Supplementary_Table_S10_TCGA_Cox_complete_with_PH.tsv")
message("Minor revision outputs written successfully.")
