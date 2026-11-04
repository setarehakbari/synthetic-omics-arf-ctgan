# Held-out real-data evaluation of synthetic augmentation

Code for the analysis added in response to Reviewer 3 (Major 2) and Reviewer 4 (points 1–3).
**Patient-level expression data are not part of this repository and must not be committed.**

## What the analysis does

For each of `N_SPLITS` stratified 75/25 splits of the 157 real subjects:

1. **Preprocessing is fitted on the training subjects only** (quantile-normalisation reference,
   log2, CoV filter 0 < CoV < 0.025) and then applied unchanged to the held-out subjects.
2. **Generators are trained on the training subjects only** (ARF, CTGAN, CTGAN with a larger update
   budget, TVAE) and each produces 200 synthetic subjects.
3. **All classifiers are scored on the same held-out real subjects.** Training sets: real-train,
   synthetic only, real-train + synthetic. Classifiers: random forest and penalised logistic
   regression (ridge, `glmnet` alpha = 0 by default; set `HOLDOUT_LR_ALPHA=0.5` for elastic net; standardisation and penalty chosen inside the training set). Metrics: accuracy,
   balanced accuracy, AUC (positive class fixed to MDD, no folding of AUC below 0.5).
4. **Dilution is recomputed inside every split**, together with two control conditions:
   Control A (real-train plus rows resampled from real-train) and Control B (same data, different
   random-forest seed).

Summaries across splits are descriptive (mean, SD, min, max). Splits share subjects, so they are
not independent replicates and the ranges are not confidence intervals.

## Requirements

* R ≥ 4.3 with `ranger`, `glmnet`, `pROC`, `dplyr`, `tibble`, `arf`
* Python ≥ 3.10 with `sdv`, `torch`, `pandas`, `numpy`

## Inputs (kept private, in `HOLDOUT_BASE_DIR`)

* `sense.filtered.cpm.Rdata` — object `sense.filtered.cpm` (genes × subjects, CPM)
* `Demographic_symptom.csv` — first column = subject ID, column `Diag` = HC / MDD
* `gene_with_label.csv` — optional, only to verify that the preprocessing reproduces the manuscript matrix

## Run order

```bash
export HOLDOUT_BASE_DIR="$HOME/synth_eval_data"   # folder that holds the inputs above

Rscript 01_make_splits.R                  # splits, seeds, train-only preprocessing (checks reproduction of the manuscript matrix)
Rscript 02a_generate_arf.R                # ARF
python3 02b_generate_ctgan_tvae.py        # CTGAN, CTGAN_long, TVAE
Rscript 03_evaluate_holdout.R             # utility + dilution + summaries
```

* `HOLDOUT_N_SPLITS` (stage 1, default 5) sets the number of splits.
* `HOLDOUT_SPLITS=1,2` (stages 2a, 3) and `--splits 1 2` (stage 2b) restrict a run to some splits, so
  several terminals can work on different splits at the same time.
* Stage 1 refuses to run if splits already exist, so that old synthetic data and results can never be
  mixed with new splits. `HOLDOUT_RESET=1 Rscript 01_make_splits.R` deletes **all** earlier outputs in
  the output folder and starts over. Do not re-run stage 1 otherwise.
* Stages 2a/2b keep existing outputs; set `HOLDOUT_OVERWRITE=1` (R) or pass `--overwrite` (Python) to redo them.
* Stage 3 re-evaluates a split when a generator file was added, replaced or changed since the last
  evaluation, and prints a WARNING if the generators were evaluated on different numbers of splits.
  Run stage 3 for the final results only after every generator has finished on every split.
* `HOLDOUT_QUICK=1` (stage 3) and `--quick` (stage 2b) are smoke tests with fewer trees/epochs.
  Their files carry a `_QUICK` / `.QUICK` suffix and never replace the real outputs. Their content is not a result.
* Every ranking comparison (observed, Control A, Control B) is repeated `HOLDOUT_CTRL_REPS` times
  (default 10) inside each split: different random-forest seeds for the observed and Control B
  rankings; a new resample and a new seed for Control A. Repetition 1 uses the seed stored in
  `seeds.csv`; the others use seeds derived from it. These repetitions measure ranking noise given the
  training data. They do not repeat the synthetic-data generation, and the controls are control
  conditions, not independent null replicates.

## Rough run time (single measurement, 2-CPU cloud container, CPU only, 117 training subjects)

| Step | Approx. time per split |
|---|---|
| ARF | 2 min |
| CTGAN (50 epochs) | 7 min |
| TVAE (300 epochs) | 9 min |
| CTGAN_long (600 epochs) | 50 min |
| Evaluation (full settings) | not yet timed |

Each CTGAN/TVAE fit has a fixed cost of about 3.5 min (data transformation of ~1,800 genes) plus about
4.6 s per CTGAN epoch and 1 s per TVAE epoch. Other machines can be slower; time one split first.
`CTGAN_long` dominates and can be run in separate terminals:
`python3 02b_generate_ctgan_tvae.py --models CTGAN_long --splits 1 2`.

## Training budgets

`generator_budget_log.csv` records epochs, batch size, updates per epoch, total updates and
wall-clock time of every fit. CTGAN uses `max(n // batch, 1)` updates per epoch and TVAE
`ceil(n / batch)`. `CTGAN_long` uses the manuscript CTGAN settings with as many updates as the
manuscript TVAE run. A CTGAN update (discriminator step plus generator step) and a TVAE update
(one joint step) are not computationally equivalent, so matched counts are a reference point, not
proof of equal training quality.

## Outputs (`reviewer_revision/R3_major2_holdout/`)

| File | Content |
|---|---|
| `seeds.csv` | collision-free seed for every (split, role) |
| `split_membership.csv`, `split_summary.csv` | which subject is train/test; sizes and class counts |
| `holdout_utility_per_split.csv`, `holdout_utility_summary.csv` | metrics on held-out real subjects |
| `holdout_utility_delta_vs_real.csv` | change versus the real-train baseline, per model and generator |
| `holdout_dilution_per_split.csv`, `holdout_dilution_summary.csv` | observed dilution and both controls |
| `holdout_dilution_within_split_spread.csv` | per split: mean, min and max of the repetitions for every comparison |
| `holdout_dilution_within_split_comparison.csv` | in how many splits the observed shift exceeds each control (means, and all repetitions) |
| `generator_budget_log.csv`, `run_info.txt` | budgets, timings, package versions |
