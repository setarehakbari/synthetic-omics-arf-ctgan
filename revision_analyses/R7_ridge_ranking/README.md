# Non-tree ranking sensitivity analysis (Reviewer 4, comment 2)

**Question answered.** The dilution and rank-agreement analyses rank genes with random-forest permutation importance.
Does the picture change when genes are ranked by a non-tree model? The SAME synthetic draws (30 per generator, plus the
five CTGAN-long draws) and the SAME controls as in `R4_repeated_draws` are re-ranked with ridge logistic regression
(genes ranked by the absolute standardized coefficient).

| Script | What it does | Time |
|---|---|---|
| `08_ridge_ranking_dilution.R` | ridge reference ranking, ridge rankings of every draw and control, metrics, summaries | about 20-40 min |

## Run
Needs `$SYNTH_EVAL_BASE_DIR/gene_with_label.csv` and the finished folder `reviewer_revision/R4_repeated_draws`
(`draw_seeds.csv` and `draws/<generator>/draw_XX.csv`). In Terminal:
```
cd revision_analyses/R7_ridge_ranking
Rscript 08_ridge_ranking_dilution.R
```
(put the script into that folder first). A first quick test (3 draws, 5 controls, results in `R7_ridge_ranking_QUICK`, not usable as results):
`RD_QUICK=1 Rscript 08_ridge_ranking_dilution.R`.

## Design
* **Ridge ranking:** logistic regression with an L2 penalty (alpha = 0), predictors standardized with the data being fitted,
  penalty chosen by stratified 5-fold cross-validation (lambda.min); genes ranked by |coefficient| (ties by gene name).
  `ridge_cv` (main) re-selects the penalty in every fit; `ridge_fixed` (sensitivity) keeps the penalty of the reference fit.
* **Reference ranking:** ridge on the 157 real subjects (CV seed 42).
* **Observed / Control A / Control B** as in the random-forest analysis. Ridge has no random component except the
  cross-validation folds, so Control B (other fold seeds) measures only the effect of the penalty choice and can be very small;
  Control A (200 resampled real rows) is the informative control here.
* **Metrics:** Spearman, Kendall, RBO (p = 0.95 and 0.98), Top-K Jaccard (25, 50, 100), mean and median absolute rank shift and
  dropout of the reference Top-50.
* Everything is descriptive and conditional on the same 157 real subjects.

## Outputs (`R7_ridge_ranking/`)
`per_ranking_all_metrics_ridge.csv`, `summary_by_condition_ridge.csv`, `table_mean_sd_range_ridge.csv`,
`exceedance_vs_controls_ridge.csv`, `generator_pairwise_comparison_ridge.csv`, `reference_ranking_ridge.csv`, `run_info_ridge.txt`.
