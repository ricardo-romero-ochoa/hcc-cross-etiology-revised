
# =============================================================================
# TCGA-LIHC lightweight module scoring and Cox modeling
#
# This script intentionally avoids GDCprepare()/SummarizedExperiment assembly,
# which can require several GB of RAM and may be killed by the operating system
# on desktop machines. It reads downloaded STAR-count files one at a time,
# extracts only ProlifHub/HepLoss genes, computes log2(CPM + 1) module scores,
# merges clinical metadata, and writes the reviewer-facing TCGA outputs.
# =============================================================================

source("R/_shared.R")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)

source("R/_shared.R")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)

  library(survival)
  library(broom)
})

ensure_dirs()
panels <- read_module_panels()
<<<<<<< HEAD
panels <- lapply(panels, function(x) toupper(unique(as.character(x))))
module_genes <- sort(unique(c(panels$ProlifHub, panels$HepLoss)))

# -----------------------------
# Generic helpers
# -----------------------------

panels <- lapply(panels, toupper)

# -----------------------------
# Helpers specific to TCGA/GDC
# -----------------------------
make_coldata_tibble <- function(se) {
  meta0 <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)

  # GDCprepare/TCGAbiolinks objects may already contain a sample_id column.
  # rownames_to_column("sample_id") fails in that case, so preserve any existing
  # field under a different name and use the SE column names as the canonical ID.
  if ("sample_id" %in% names(meta0)) {
    names(meta0)[names(meta0) == "sample_id"] <- "sample_id_coldata"
  }
  names(meta0) <- make.unique(names(meta0), sep = "_")

  meta <- tibble::rownames_to_column(meta0, var = "sample_id") |>
    tibble::as_tibble()

  # Some SE objects use simple rownames but expression colnames carry the TCGA barcode.
  # If possible, force metadata IDs to match the assay colnames exactly.
  if (nrow(meta) == ncol(se)) {
    meta$sample_id <- colnames(se)
  }
  meta
}

>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

<<<<<<< HEAD
coalesce_cols <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep(NA, nrow(df)))
  out <- df[[cols[1]]]
  if (length(cols) > 1) {
    for (cc in cols[-1]) out <- dplyr::coalesce(out, df[[cc]])
  }
  out
}

normalize_barcode <- function(x) {
  x <- as.character(x)
  x <- stringr::str_split_fixed(x, ";", 2)[, 1]
  x <- stringr::str_trim(x)
  x
}

infer_tcga_tissue_from_barcode <- function(sample_id, sample_type_text = NA_character_) {
  sample_id <- as.character(sample_id)
  code <- suppressWarnings(substr(sample_id, 14, 15))
  by_code <- dplyr::case_when(
    code %in% c("01", "02", "03", "05", "06", "07", "08", "09") ~ "tumor",
    code %in% c("10", "11", "12", "13", "14") ~ "non_tumor",
    TRUE ~ NA_character_
  )

  txt <- as.character(sample_type_text)
  by_text <- dplyr::case_when(
    stringr::str_detect(txt, stringr::regex("solid tissue normal|normal", TRUE)) ~ "non_tumor",
    stringr::str_detect(txt, stringr::regex("primary tumor|tumou?r|carcinoma|cancer", TRUE)) ~ "tumor",
    TRUE ~ NA_character_
  )
  dplyr::coalesce(by_code, by_text)
}

locate_gdc_files <- function(res, root = "GDCdata/TCGA-LIHC") {
  all_files <- list.files(root, pattern = "star_gene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) {
    stop("No downloaded STAR-count files found under ", root,
         ". Run GDCdownload(query) first or check the GDCdata directory.")
  }

  file_col <- first_existing_col(res, c("file_name", "filename", "file.name"))
  id_col <- first_existing_col(res, c("file_id", "id", "file_id.x"))
  if (is.na(file_col)) stop("Could not find file_name column in TCGAbiolinks getResults(query).")

  res$file_path <- all_files[match(res[[file_col]], basename(all_files))]

  if (anyNA(res$file_path) && !is.na(id_col)) {
    miss <- which(is.na(res$file_path))
    for (i in miss) {
      fid <- as.character(res[[id_col]][i])
      hit <- all_files[stringr::str_detect(all_files, fixed(fid))]
      if (length(hit) > 0) res$file_path[i] <- hit[[1]]
    }
  }

  missing_n <- sum(is.na(res$file_path))
  if (missing_n > 0) {
    warning("Could not locate ", missing_n, " downloaded TCGA count files. These samples will be skipped.")
  }
  res |> dplyr::filter(!is.na(.data$file_path))
}

read_one_star_file <- function(path, sample_id, module_genes) {
  header <- names(data.table::fread(path, nrows = 0, showProgress = FALSE))
  gene_id_col <- first_existing_col(as.data.frame(matrix(ncol = length(header), nrow = 0, dimnames = list(NULL, header))),
                                    c("gene_id", "Geneid", "gene", "ensembl_gene_id"))
  gene_name_col <- first_existing_col(as.data.frame(matrix(ncol = length(header), nrow = 0, dimnames = list(NULL, header))),
                                      c("gene_name", "gene_symbol", "external_gene_name", "hgnc_symbol"))
  count_col <- first_existing_col(as.data.frame(matrix(ncol = length(header), nrow = 0, dimnames = list(NULL, header))),
                                  c("unstranded", "expected_count", "counts"))
  tpm_col <- first_existing_col(as.data.frame(matrix(ncol = length(header), nrow = 0, dimnames = list(NULL, header))),
                                c("tpm_unstranded", "TPM", "tpm"))

  if (is.na(gene_name_col)) {
    stop("No gene_name/gene_symbol column found in ", path, ". Columns: ", paste(header, collapse = ", "))
  }
  if (is.na(count_col) && is.na(tpm_col)) {
    stop("No count or TPM column found in ", path, ". Columns: ", paste(header, collapse = ", "))
  }

  value_col <- if (!is.na(count_col)) count_col else tpm_col
  select_cols <- unique(stats::na.omit(c(gene_id_col, gene_name_col, value_col)))
  dat <- data.table::fread(path, select = select_cols, showProgress = FALSE) |>
    tibble::as_tibble()

  names(dat)[names(dat) == gene_name_col] <- "gene_name"
  names(dat)[names(dat) == value_col] <- "value"
  if (!is.na(gene_id_col) && gene_id_col %in% names(dat)) names(dat)[names(dat) == gene_id_col] <- "gene_id"
  if (!"gene_id" %in% names(dat)) dat$gene_id <- NA_character_

  dat <- dat |>
    dplyr::mutate(
      gene = toupper(as.character(.data$gene_name)),
      value = suppressWarnings(as.numeric(.data$value)),
      is_summary_row = stringr::str_detect(as.character(.data$gene_id), "^N_")
    )

  using_counts <- !is.na(count_col)
  lib_size <- if (using_counts) {
    sum(dat$value[!dat$is_summary_row & !is.na(dat$value)], na.rm = TRUE)
  } else {
    NA_real_
  }

  mod <- dat |>
    dplyr::filter(!.data$is_summary_row, .data$gene %in% module_genes) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(raw_value = sum(.data$value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      sample_id = sample_id,
      library_size = lib_size,
      value_type = ifelse(using_counts, "counts", "tpm")
    )

  mod
=======
infer_tcga_tissue <- function(meta) {
  # Prefer the TCGA barcode sample code because it is stable across TCGAbiolinks versions:
  # 01 = primary tumor; 11 = solid tissue normal. Fall back to text metadata.
  sample_code <- suppressWarnings(substr(meta$sample_id, 14, 15))
  tissue_by_code <- dplyr::case_when(
    sample_code %in% c("01", "02", "03", "05", "06", "07", "08", "09") ~ "tumor",
    sample_code %in% c("10", "11", "12", "13", "14") ~ "non_tumor",
    TRUE ~ NA_character_
  )

  text_cols <- intersect(
    c("sample_type", "definition", "shortLetterCode", "tissue_type", "sample", "sample_id_coldata"),
    names(meta)
  )
  text <- rep("", nrow(meta))
  if (length(text_cols) > 0) {
    text <- apply(meta[, text_cols, drop = FALSE], 1, function(x) {
      paste(as.character(x[!is.na(x)]), collapse = " | ")
    })
  }

  tissue_by_text <- dplyr::case_when(
    stringr::str_detect(text, stringr::regex("solid tissue normal|normal tissue|adjacent|non[- ]?tumou?r|non[- ]?cancer", TRUE)) ~ "non_tumor",
    stringr::str_detect(text, stringr::regex("primary tumor|tumou?r|carcinoma|cancer", TRUE)) ~ "tumor",
    TRUE ~ NA_character_
  )

  dplyr::coalesce(tissue_by_code, tissue_by_text)
}

collapse_counts_to_symbols <- function(expr_counts, rowdata) {
  rd <- as.data.frame(rowdata, stringsAsFactors = FALSE)
  sym_col <- first_existing_col(rd, c("gene_name", "external_gene_name", "gene_symbol", "hgnc_symbol"))
  if (is.na(sym_col)) {
    stop("Could not find a gene-symbol column in TCGA rowData. Columns found: ", paste(names(rd), collapse = ", "))
  }

  sym <- toupper(as.character(rd[[sym_col]]))
  ok <- !is.na(sym) & sym != "" & !stringr::str_detect(sym, "^ENSG")
  expr_counts <- expr_counts[ok, , drop = FALSE]
  sym <- sym[ok]
  storage.mode(expr_counts) <- "numeric"
  rowsum(expr_counts, group = sym, reorder = FALSE)
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
}

clean_clinical_for_merge <- function(clin) {
  clin <- tibble::as_tibble(clin)
  if (!"submitter_id" %in% names(clin)) {
    stop("TCGA clinical table lacks submitter_id; cannot merge clinical metadata.")
  }
  clin |>
    dplyr::mutate(patient = .data$submitter_id) |>
    dplyr::distinct(.data$patient, .keep_all = TRUE)
}

<<<<<<< HEAD
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

=======
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
fit_lm_associations <- function(scores2) {
  candidates <- c("ajcc_pathologic_stage", "tumor_grade", "grade", "gender", "sex", "age_at_index")
  score_names <- c("ProlifHubScore", "HepLossScore", "HCCStateScore")
  assoc <- list()

  for (v in intersect(candidates, names(scores2))) {
    for (s in score_names) {
      dat <- scores2 |>
        dplyr::filter(.data$tissue == "tumor", !is.na(.data[[v]]), !is.na(.data[[s]]))
      if (nrow(dat) < 10) next
      if (!is.numeric(dat[[v]])) {
        dat[[v]] <- as.factor(dat[[v]])
        if (nlevels(dat[[v]]) < 2) next
      }
      fit <- tryCatch(stats::lm(stats::reformulate(v, response = s), data = dat), error = function(e) NULL)
      if (!is.null(fit)) {
        assoc[[paste(v, s, sep = "__")]] <- broom::tidy(fit) |>
          dplyr::mutate(variable = v, score = s, n = nrow(dat))
      }
    }
  }

<<<<<<< HEAD
  if (length(assoc) == 0) return(tibble::tibble())
=======
  if (length(assoc) == 0) {
    return(tibble::tibble())
  }
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
  dplyr::bind_rows(assoc)
}

fit_survival_models <- function(scores2) {
<<<<<<< HEAD
=======
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

>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
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

  readr::write_tsv(survdat, "results/tables/tcga_lihc_survival_model_input.tsv")

  if (nrow(survdat) < 30 || length(unique(survdat$event)) < 2) {
    message("[06] Skipping Cox models: insufficient survival events or samples.")
    return(tibble::tibble())
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
        message("[06] Cox model failed for ", score_name, " / ", model_name, ": ", conditionMessage(e))
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
<<<<<<< HEAD
  dplyr::bind_rows(
=======
  out <- dplyr::bind_rows(
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
    lapply(scores_to_test, fit_one, covars = character(), model_name = "score_only"),
    lapply(scores_to_test, fit_one, covars = cov_age_sex, model_name = "age_sex_adjusted"),
    lapply(scores_to_test, fit_one, covars = cov_age_sex_stage, model_name = "age_sex_stage_adjusted")
  )
<<<<<<< HEAD
}

# -----------------------------
# Query/download metadata only; read expression one file at a time
# -----------------------------
message("[06] Querying TCGA-LIHC STAR-count files.")
=======
  out
}

# -----------------------------
# Download and prepare TCGA-LIHC
# -----------------------------
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
query <- TCGAbiolinks::GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
TCGAbiolinks::GDCdownload(query)
<<<<<<< HEAD

res <- TCGAbiolinks::getResults(query) |>
  tibble::as_tibble()

case_col <- first_existing_col(res, c(
  "cases", "case_submitter_id", "cases.submitter_id", "submitter_id",
  "sample_submitter_id", "sample.submitter_id", "cases.samples.submitter_id"
))
sample_type_col <- first_existing_col(res, c(
  "sample_type", "cases.samples.sample_type", "samples.sample_type", "sample_type_id"
))
file_col <- first_existing_col(res, c("file_name", "filename", "file.name"))
file_id_col <- first_existing_col(res, c("file_id", "id", "file_id.x"))

if (is.na(case_col)) {
  stop("Could not find a TCGA barcode/case column in getResults(query). Columns found: ", paste(names(res), collapse = ", "))
}
if (is.na(file_col)) {
  stop("Could not find file_name in getResults(query). Columns found: ", paste(names(res), collapse = ", "))
}

res <- res |>
  dplyr::mutate(
    case_barcode_raw = normalize_barcode(.data[[case_col]]),
    file_identifier = if (!is.na(file_id_col)) as.character(.data[[file_id_col]]) else as.character(.data[[file_col]]),
    patient = substr(.data$case_barcode_raw, 1, 12),
    sample_id_base = dplyr::if_else(
      nchar(.data$case_barcode_raw) >= 16,
      .data$case_barcode_raw,
      paste0(.data$case_barcode_raw, "_", substr(.data$file_identifier, 1, 8))
    ),
    sample_id = make.unique(.data$sample_id_base, sep = "_dup"),
    sample_type_text = if (!is.na(sample_type_col)) as.character(.data[[sample_type_col]]) else NA_character_,
    tissue = infer_tcga_tissue_from_barcode(.data$case_barcode_raw, .data$sample_type_text)
  ) |>
  locate_gdc_files()

readr::write_tsv(res, "results/tables/tcga_lihc_gdc_file_manifest.tsv")

message("[06] Reading ", nrow(res), " STAR-count files one at a time; extracting ", length(module_genes), " module genes.")
expr_long <- purrr::map2_dfr(
  res$file_path,
  res$sample_id,
  function(path, sid) {
    tryCatch(
      read_one_star_file(path, sid, module_genes),
      error = function(e) {
        warning("Failed to read ", path, " for ", sid, ": ", conditionMessage(e))
        tibble::tibble(gene = character(), raw_value = numeric(), sample_id = character(), library_size = numeric(), value_type = character())
      }
    )
  }
)

if (nrow(expr_long) == 0) {
  stop("No module-gene expression values were extracted from TCGA-LIHC files.")
}

readr::write_tsv(expr_long, "results/tables/tcga_lihc_module_gene_raw_values.tsv")

# Complete absent module genes as zero counts/TPM before scoring.
expr_long2 <- tidyr::expand_grid(sample_id = unique(res$sample_id), gene = module_genes) |>
  dplyr::left_join(expr_long, by = c("sample_id", "gene")) |>
  dplyr::left_join(res |> dplyr::select(.data$sample_id, .data$file_path), by = "sample_id") |>
  dplyr::group_by(.data$sample_id) |>
  dplyr::mutate(
    library_size = dplyr::coalesce(.data$library_size, .data$library_size[which(!is.na(.data$library_size))[1]]),
    value_type = dplyr::coalesce(.data$value_type, .data$value_type[which(!is.na(.data$value_type))[1]], "counts"),
    raw_value = dplyr::coalesce(.data$raw_value, 0)
  ) |>
  dplyr::ungroup()

expr_df <- expr_long2 |>
  dplyr::mutate(
    expression = dplyr::if_else(
      .data$value_type == "counts",
      log2((.data$raw_value / .data$library_size) * 1e6 + 1),
      log2(.data$raw_value + 1)
    )
  ) |>
  dplyr::select(.data$gene, .data$sample_id, .data$expression) |>
  tidyr::pivot_wider(names_from = .data$sample_id, values_from = .data$expression, values_fill = 0)

expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- expr_df$gene
storage.mode(expr_mat) <- "numeric"

meta <- res |>
  dplyr::select(.data$sample_id, .data$patient, .data$tissue, .data$sample_type_text, .data$file_path) |>
  dplyr::distinct(.data$sample_id, .keep_all = TRUE) |>
  dplyr::arrange(match(.data$sample_id, colnames(expr_mat)))
expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]
=======
se <- TCGAbiolinks::GDCprepare(query)

expr_counts <- SummarizedExperiment::assay(se)
expr_counts <- collapse_counts_to_symbols(expr_counts, SummarizedExperiment::rowData(se))
expr_log <- log2(edgeR::cpm(expr_counts, log = FALSE) + 1)
rownames(expr_log) <- toupper(rownames(expr_log))

meta <- make_coldata_tibble(se) |>
  dplyr::mutate(tissue = infer_tcga_tissue(dplyr::cur_data_all()))

# Keep only samples that can be interpreted as tumor or normal/non-tumor.
meta <- meta |>
  dplyr::filter(.data$sample_id %in% colnames(expr_log)) |>
  dplyr::arrange(match(.data$sample_id, colnames(expr_log)))
expr_log <- expr_log[, meta$sample_id, drop = FALSE]
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4

label_audit <- meta |>
  dplyr::count(.data$tissue, name = "n") |>
  dplyr::mutate(dataset = "TCGA-LIHC", .before = 1)
readr::write_tsv(label_audit, "results/tables/tcga_lihc_tissue_label_audit.tsv")

if (!all(c("tumor", "non_tumor") %in% unique(meta$tissue))) {
  message("[06] Warning: TCGA-LIHC did not contain both tumor and non_tumor labels after parsing. Continuing with available labels.")
}

<<<<<<< HEAD
scores <- score_modules(expr_mat, panels$ProlifHub, panels$HepLoss) |>
  dplyr::left_join(meta, by = "sample_id")
readr::write_tsv(scores, "results/tables/tcga_lihc_module_scores.tsv")

# -----------------------------
# Clinical metadata and models
# -----------------------------
message("[06] Downloading/reading TCGA-LIHC clinical metadata.")
=======
scores <- score_modules(expr_log, panels$ProlifHub, panels$HepLoss) |>
  dplyr::left_join(meta, by = "sample_id")
readr::write_tsv(scores, "results/tables/tcga_lihc_module_scores.tsv")

# Clinical metadata
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
clin <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LIHC", type = "clinical") |>
  tibble::as_tibble()
readr::write_tsv(clin, "results/tables/tcga_lihc_clinical_raw.tsv")
clin2 <- clean_clinical_for_merge(clin)

scores2 <- scores |>
  dplyr::mutate(patient = substr(.data$sample_id, 1, 12)) |>
  dplyr::left_join(clin2, by = "patient")
readr::write_tsv(scores2, "results/tables/tcga_lihc_module_scores_with_clinical.tsv")

assoc <- fit_lm_associations(scores2)
readr::write_tsv(assoc, "results/tables/tcga_lihc_clinicopathologic_associations.tsv")

cox <- fit_survival_models(scores2)
readr::write_tsv(cox, "results/tables/tcga_lihc_survival_cox_models.tsv")

<<<<<<< HEAD
message("[06] TCGA-LIHC lightweight pipeline completed. Outputs written to results/tables/.")
=======
message("[06] TCGA-LIHC pipeline completed. Outputs written to results/tables/.")
>>>>>>> de0c748d3b558ba656b9e43a99b5bec165a230a4
