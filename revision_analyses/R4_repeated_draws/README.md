# Repeated synthetic-data draws (Reviewer 4 point 1, Reviewer 3 Major 1)

**Question answered.** The manuscript reports one synthetic data set per generator. Here the complete
generation procedure is repeated 30 times per generator on the full real data (157 subjects x 1,783 genes),
each time with a new seed, a new model fit and a new sample of 200 synthetic subjects. Every draw is then
scored with the same ranking pipeline as Table 2 and compared with two real-data-only controls.

| Step | Script | What it does | Expected time (estimate) |
|---|---|---|---|
| 0 | `05_prepare_repeated_draws.R` | checks the real data, writes `draw_seeds.csv` (never overwritten) | seconds |
| 1 | `05a_draws_arf.R` | 30 ARF draws | 25-40 min |
| 2 | `05b_draws_ctgan_tvae.py` | 30 CTGAN + 30 TVAE draws (draw by draw, resumable) | CTGAN ~2-3 h, TVAE ~6-8 h |
| 3 | `05c_evaluate_repeated_draws.R` | metrics, controls, summaries, figure | 15-25 min |

Time the first draw and multiply by 30. Every script skips finished draws, so it can be stopped and restarted.

## Run (Mac; keep the Mac plugged in, awake, and Terminal open)
Put the five files in one folder. R scripts: in R/RStudio run `source(file.choose())` and pick the file (or `Rscript <file>` in Terminal).
```
05_prepare_repeated_draws.R        # R: seconds
05a_draws_arf.R                    # R: ~30 min
```
Python script: in Terminal type `cd ` (with a space), drag the folder into the window, press Enter, then
```
caffeinate -i python3 05b_draws_ctgan_tvae.py
```
When it has finished (the last line says "CTGAN/TVAE finished."), run `05c_evaluate_repeated_draws.R` in R.
The scripts find the data through `$SYNTH_EVAL_BASE_DIR` (or the environment variable `RD_BASE_DIR`).

Optional, slow: `python3 05b_draws_ctgan_tvae.py --first 1 --last 5 --models CTGAN_long`
(CTGAN with as many updates as TVAE: 900 updates at n = 157; roughly 30-35 min per draw).
Smoke tests (separate folders, never usable as results): `RD_QUICK=1` for R, `--quick` for Python.

## Design
* **One draw = full refit.** ARF: `adversarial_rf(100 trees)` -> `forde(alpha 0.1, epsilon 0.1)` -> `forge(200)`.
  CTGAN: 50 epochs, batch 256, PAC 1 (1 update/epoch at n = 157, i.e. 50 updates).
  TVAE: 300 epochs, batch 64, embedding 128, (128,128), L2 1e-5, loss factor 2 (3 updates/epoch, 900 updates).
* **Reference ranking:** random forest on the real data, permutation importance, 3000 trees, seed 42 (as in Table 2).
* **Observed:** same random forest on real + synthetic draw (own RF seed per draw).
* **Control A:** real + 200 rows resampled with replacement from the real data, same HC/MDD counts as that draw (no generator).
* **Control B:** real data only, re-ranked with 200 other RF seeds (pure random-forest noise).
* **Metrics:** Spearman, Kendall, RBO (p = 0.95; p = 0.98 as sensitivity), Top-K Jaccard (K = 25/50/100),
  mean/median absolute rank shift and % dropout of the real Top-50, mean/median signed shift.
* The manuscript's single draws (`synthetic_rna_arf.csv`, `synthetic_ctgan.csv`, `synthetic_rna_tvae.csv`) are
  re-scored as "draw 0": the console says whether Table 2 is reproduced (mean absolute shift 418.36 / 481.66 / 413.38).

## Outputs (`results/`)
`table_mean_sd_range.csv` (mean (SD) [min to max] per condition: the table for the manuscript),
`summary_by_condition.csv` (mean, SD, min, 2.5/50/97.5 percentiles, max), `exceedance_vs_controls.csv`
(draws beyond every control value; empirical p = (1 + #controls as or more extreme)/(B + 1)),
`generator_pairwise_comparison.csv` (mean difference with bootstrap 95% CI, probability of superiority, Wilcoxon p),
`per_draw_all_metrics.csv`, `synthetic_class_balance.csv`, `generation_log_all.csv` (budgets, run times, versions),
`fig_repeated_draws.png/.pdf`, `run_info.txt`.

## Interpretation limits (state them in the manuscript)
* All draws are conditional on the same 157 real subjects. The intervals show run-to-run variation of a generator and
  its separation from ranking noise; they are not confidence intervals for the real population.
* Control A repeats existing individuals (re-weighting), Control B changes no data; neither is a new independent sample.
* Pairwise comparisons are descriptive and not adjusted for multiplicity.

## Privacy
`draws/` and `results/per_draw_all_metrics.csv` derive from patient data: keep `draws/`, `*.csv` with gene profiles and
`gene_with_label.csv` out of any public repository (code only).
