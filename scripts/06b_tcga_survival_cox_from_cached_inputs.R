# =============================================================================
# Rebuild TCGA-LIHC Cox survival table from cached module-score/clinical table
#
# Use when scripts/06_tcga_lihc_pipeline.R already produced
# results/tables/tcga_lihc_module_scores_with_clinical.tsv but the final Cox table
# needs to be regenerated without re-downloading TCGA/GDC data.
# =============================================================================

source("R/_shared.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(broom)
})

ensure_dirs()
input_file <- "results/tables/tcga_lihc_module_scores_with_clinical.tsv"
output_file <- "results/tables/tcga_lihc_survival_cox_models.tsv"
model_input_file <- "results/tables/tcga_lihc_survival_model_input.tsv"

if (!file.exists(input_file)) {
  stop("Missing ", input_file, ". Run scripts/06_tcga_lihc_pipeline.R first.")
}

scores2 <- readr::read_tsv(input_file, show_col_types = FALSE)

coalesce_cols <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep(NA, nrow(df)))
  out <- df[[cols[1]]]
  if (length(cols) > 1) {
    for (cc in cols[-1]) out <- dplyr::coalesce(out, df[[cc]])
  }
  out
}

collapse_stage <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "STAGE I[^V]|STAGE IA|STAGE IB|^I$|^IA$|^IB$") ~ "I",
    stringr::str_detect(x, "STAGE II|^II$|^IIA$|^IIB$") ~ "II",
    stringr::str_detect(x, "STAGE III|^III$|^IIIA$|^IIIB$|^IIIC$") ~ "III",
    stringr::str_detect(x, "STAGE IV|^IV$|^IVA$|^IVB$") ~ "IV",
    TRUE ~ NA_character_
  )
}

survdat <- scores2 |>
  dplyr::filter(.data$tissue == "tumor") |>
  dplyr::mutate(
    vital_status_tmp = as.character(coalesce_cols(dplyr::cur_data_all(), c("vital_status.y", "vital_status.x", "vital_status"))),
    days_to_death_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("days_to_death.y", "days_to_death.x", "days_to_death")))),
    days_to_follow_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("days_to_last_follow_up.y", "days_to_last_follow_up.x", "days_to_last_follow_up")))),
    age_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("age_at_index.y", "age_at_index.x", "age_at_index")))),
    gender_tmp = as.factor(coalesce_cols(dplyr::cur_data_all(), c("gender.y", "gender.x", "gender", "sex"))),
    stage_raw = as.character(coalesce_cols(dplyr::cur_data_all(), c("ajcc_pathologic_stage.y", "ajcc_pathologic_stage.x", "ajcc_pathologic_stage"))),
    stage_collapsed = factor(collapse_stage(stage_raw), levels = c("I", "II", "III", "IV")),
    time_days = ifelse(!is.na(days_to_death_tmp), days_to_death_tmp, days_to_follow_tmp),
    event = ifelse(stringr::str_detect(vital_status_tmp, stringr::regex("dead", ignore_case = TRUE)), 1, 0),
    time_months = time_days / 30.44
  ) |>
  dplyr::filter(is.finite(.data$time_months), .data$time_months > 0)

readr::write_tsv(survdat, model_input_file)

if (nrow(survdat) < 30 || length(unique(survdat$event)) < 2) {
  stop("Insufficient TCGA-LIHC survival events or samples after filtering; cannot fit Cox models.")
}

has_age <- sum(!is.na(survdat$age_tmp)) >= 100
has_gender <- nlevels(droplevels(survdat$gender_tmp)) > 1
has_stage <- nlevels(droplevels(survdat$stage_collapsed)) > 1

cov_age_sex <- character()
if (has_age) cov_age_sex <- c(cov_age_sex, "age_tmp")
if (has_gender) cov_age_sex <- c(cov_age_sex, "gender_tmp")

cov_age_sex_stage <- cov_age_sex
if (has_stage) cov_age_sex_stage <- c(cov_age_sex_stage, "stage_collapsed")

fit_one <- function(score_name, covars = character(), model_name = "score_only") {
  rhs <- paste0("scale(", score_name, ")")
  if (length(covars) > 0) rhs <- paste(rhs, paste(covars, collapse = " + "), sep = " + ")
  f <- stats::as.formula(paste("survival::Surv(time_months, event) ~", rhs))
  fit <- tryCatch(
    survival::coxph(f, data = survdat),
    error = function(e) {
      message("Cox model failed for ", score_name, " / ", model_name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    dplyr::mutate(
      model = model_name,
      score = score_name,
      formula = paste(deparse(f), collapse = " "),
      n = fit$n,
      events = fit$nevent,
      .before = 1
    )
}

scores_to_test <- c("ProlifHubScore", "HepLossScore", "HCCStateScore")
cox <- dplyr::bind_rows(
  lapply(scores_to_test, fit_one, covars = character(), model_name = "score_only"),
  lapply(scores_to_test, fit_one, covars = cov_age_sex, model_name = "age_sex_adjusted"),
  lapply(scores_to_test, fit_one, covars = cov_age_sex_stage, model_name = "age_sex_stage_adjusted")
)

if (nrow(cox) == 0) {
  stop("All TCGA-LIHC Cox models failed; inspect ", model_input_file)
}

readr::write_tsv(cox, output_file)
message("Wrote ", output_file)
