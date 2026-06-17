# TCGA-LIHC memory-safe pipeline

If `scripts/06_tcga_lihc_pipeline.R` is killed immediately after `GDCdownload()` reports that all files are present, the failure is almost certainly an operating-system memory kill during full TCGA object assembly. In that situation R may terminate without writing `results/logs/ERROR_06_tcga_lihc_pipeline.R.log`, because the process is killed externally rather than stopped by an R exception.

The default `scripts/06_tcga_lihc_pipeline.R` in this revision avoids this by not calling `GDCprepare()`. It reads each downloaded `star_gene_counts.tsv` file one at a time, extracts only the ProlifHub/HepLoss genes needed for module scoring, computes log2(CPM + 1), merges clinical metadata, and writes:

```text
results/tables/tcga_lihc_module_scores.tsv
results/tables/tcga_lihc_module_scores_with_clinical.tsv
results/tables/tcga_lihc_survival_model_input.tsv
results/tables/tcga_lihc_survival_cox_models.tsv
```

Run:

```bash
Rscript scripts/06_tcga_lihc_pipeline.R
```

Then regenerate downstream manuscript objects:

```bash
Rscript scripts/11_make_revision_figures.R
Rscript scripts/12_make_manuscript_tables.R
Rscript scripts/13_summarize_revision_results.R
Rscript scripts/14_make_supplementary_figures.R
Rscript scripts/15_make_supplementary_tables.R
Rscript scripts/16_revision_repository_audit.R --strict
```

The legacy full SummarizedExperiment implementation is preserved as:

```text
scripts/06_tcga_lihc_pipeline_full_se_legacy.R
```

Use the legacy script only on a machine with sufficient RAM.
