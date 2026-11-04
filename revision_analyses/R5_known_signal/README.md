# Known-signal benchmark for the dilution metric (Reviewer 4, Comment 4.1)

**Question.** In data where the informative genes are known, (i) does augmentation with ARF / CTGAN / TVAE
remove true signal from the Top-50, (ii) does the dilution of the reference Top-50 track that loss, and
(iii) does it add information beyond rank-agreement and overlap measures (Spearman, RBO, Top-K Jaccard)?

**Data.** 10 simulated data sets from `npdro::createSimulation2` (B. McKinney's settings), 157 subjects
(79 / 78 per class) x 1,800 genes, 5% signal genes = 45 main-effect + 45 interaction-only genes
(main.bias = interaction.bias = 1.5, hi.cor = 0.95, lo.cor = 0.2, "main-interactionScalefree", pct.mixed = 0.5).
The simulated data contain no patient information: this whole folder (data, draws, results) can be made public.

**Pipeline** (same generators, settings, ranking and metrics as the main analysis)
| Step | Script | What it does | Time on a Mac (estimate) |
|---|---|---|---|
| 0 | `06_make_simulations.R` | already done (sims/, truth.csv, seeds.csv); re-run only to regenerate (needs `remotes::install_github("insilico/npdro")`) | ~10 min |
| 1 | `06a_sim_arf.R` | one ARF data set (200 subjects) per simulation | ~10 min |
| 2 | `06b_sim_ctgan_tvae.py` | one CTGAN and one TVAE data set per simulation | CTGAN ~1 h, TVAE ~2-3 h |
| 3 | `06c_sim_evaluate.R` | rankings, controls, truth-based metrics, tables, figures | ~20-30 min |
| 4 | `06d_sim_sensitivity_shift_only.R` | sensitivity analysis: incremental R² with only the mean and median absolute rank shift as dilution summaries (Top-50 dropout is an exact function of Top-50 Jaccard similarity) | ~2 min |
| 5 | `06e_sim_sensitivity_group_weights.R` | sensitivity analysis: the incremental R² comparison with equal total weight for every condition x generator within a simulation (Control B has 30 rankings per simulation, so otherwise it dominates the rows) | ~4 min |

Put this folder at `the R5_known_signal folder` (or set `KS_DIR`). In Terminal:
```
cd revision_analyses/R5_known_signal
Rscript 06a_sim_arf.R
caffeinate -i python3 06b_sim_ctgan_tvae.py
Rscript 06c_sim_evaluate.R
```
Finished sims are skipped, so every script can be stopped and restarted. Smoke tests: `KS_QUICK=1` (R) or `--quick` (Python).

**Conditions per simulation:** Reference (RF permutation importance, 3000 trees); Observed (simulation + 200 synthetic subjects);
Control A (simulation + 200 resampled rows, same class counts); Control B (30 other RF seeds); Control C (simulation + 200
resampled rows with PERMUTED labels = label-permuted positive control, 3 replicates; the 157 original rows are unchanged).

**Truth-based metrics:** true main-effect genes in the Top-50 (and loss relative to the reference), true-signal dropout, noise genes
entering the Top-50, AUROC of the importance scores for main-effect vs noise genes. Interaction-only genes are reported separately
(random-forest permutation importance rarely recovers them).

**Outputs (`results/`):** `per_ranking_metrics.csv`, `table_condition_summary.csv`, `table_paired_vs_controls.csv`,
`table_metric_vs_truth_loss.csv` (within-simulation Spearman between each summary and the loss of true signal),
`table_incremental_R2.csv` (does dilution add explanatory power beyond Spearman + RBO + Jaccard, and vice versa),
`table_detection_auc.csv`, `fig_truth_recovery.png`, `fig_loss_vs_dilution.png`, `run_info.txt`.

**Limits.** 10 simulated data sets and one synthetic data set per generator and simulation; all summaries are descriptive
(means over simulations with bootstrap intervals over simulations). The simulation is one scenario (n = 157, 5% signals); real
expression data differ. Results show whether dilution tracks known-signal loss in this setting, not that it does so in general.
