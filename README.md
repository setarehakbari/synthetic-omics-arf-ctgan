# Multidimensional evaluation of synthetic gene-expression data with a dilution metric

Code for the manuscript *"Multidimensional Evaluation of Synthetic Gene Expression Data with a Dilution Metric for
Feature Ranking Preservation"* (BMC Bioinformatics, under revision).

**This repository contains code only.** The patient-level RNA-seq data and every synthetic data set derived from them are
not shared (data-use restrictions; see the *Availability of data and materials* section of the manuscript).
The known-signal benchmark (`revision_analyses/R5_known_signal`) uses simulated data that contain no patient information
and can be regenerated with the scripts.

The study evaluates synthetic gene-expression data generated with the Adversarial Random Forest (ARF), the Conditional
Tabular GAN (CTGAN) and the Tabular Variational Autoencoder (TVAE). The evaluation framework considers six complementary
dimensions: predictive utility, statistical fidelity, privacy-related distinguishability, feature-importance alignment,
dilution of real-data feature rankings, and bootstrap feature-selection stability. In the revision, dilution is compared
with real-data-only controls and across repeated synthetic draws, utility is also evaluated on held-out real subjects, and
a benchmark with known informative genes is used.

## Repository structure

| Folder | Content |
|---|---|
| `main_pipeline/` | main analysis (see below) |
| `generators/` | synthetic-data generation: `ARF_generate.R` (arf package), `CTGAN_generate.py` and `TVAE_generate.py` (SDV), with the settings and seeds stated in the manuscript |
| `figures/` | scripts that redraw the manuscript figures and the workflow diagram |
| `revision_analyses/R3_dilution_controls_50rep` | first set of real-data-only controls for dilution (50 replicates) |
| `revision_analyses/R3_heldout_utility` | held-out real-subject utility analysis (preprocessing and generators fitted on training subjects only) |
| `revision_analyses/R4_repeated_draws` | 30 independent synthetic draws per generator, matched resampling controls (Control A) and RF-seed controls (Control B) |
| `revision_analyses/R5_known_signal` | known-signal benchmark in 10 simulated data sets, sensitivity analyses |
| `revision_analyses/R6_privacy_fold_auc` | domain-separator AUC within each validation fold (random forest; ridge logistic regression) |
| `revision_analyses/R7_ridge_ranking` | repeated-draw dilution analysis with rankings from ridge logistic regression |
| `quantum/` | exploratory quantum–adversarial robustness analysis (Additional file 1); real data only, expects the real data in a `data/` folder that is not part of this repository |
| `tools/check_before_push.sh` | checks that no data-like file or local path is about to be committed |

Each `revision_analyses/*` folder has its own README (where present) with the run order.

### `main_pipeline/`

* `preprocess_real_data.R`: reproduces quantile normalization, log2 transformation, coefficient-of-variation (CoV)
  filtering, zero-variance feature removal, and construction of the final sample-by-gene analysis matrix, for researchers
  with authorized access to the source data.
* `evaluate_all_generators.R`: primary evaluation pipeline for ARF, CTGAN and TVAE (utility, fidelity, privacy-related
  distinguishability, RF and NPDR feature alignment, dilution, bootstrap stability).
* `Table2_rank_metrics_with_TVAE.R`: ranking-agreement and dilution measures for the three generators. The single-dataset
  ranking-agreement and dilution values (Table 3 of the manuscript) come from this analysis.
* `topk_sensitivity_analysis.R`: evaluates the sensitivity of feature-alignment and dilution measures to K = 25, 50 and 100
  using the fixed full feature rankings from the primary analysis.

## Feature-ranking and dilution analysis

The feature-ranking analysis compares feature-importance rankings derived from the real data with rankings obtained from
synthetic and real-plus-synthetic datasets. The following complementary measures are evaluated:

* Spearman's rank correlation
* Kendall's rank correlation
* Rank-biased overlap (RBO; p = 0.95, with p = 0.98 as a sensitivity analysis)
* Top-K Jaccard similarity
* Mean absolute rank shift
* Median absolute rank shift
* Top-K dropout

The primary manuscript analysis uses K = 50. Mean absolute rank shift is used as the primary summary of dilution magnitude,
while signed rank displacement can be retained as a secondary directional summary. Dilution is calculated using random
forest permutation-importance rankings; NPDR is evaluated separately as an alternative feature-relevance method, and ridge
logistic regression is used in a ranking-method sensitivity analysis (`revision_analyses/R7_ridge_ranking`).

### Top-K sensitivity analysis

Feature-alignment and dilution analyses were additionally evaluated at K = 25 and K = 100, with K = 50 retained as the
primary setting. `main_pipeline/topk_sensitivity_analysis.R` reuses the fixed full RF and NPDR feature rankings generated in
the primary analysis; therefore, no RF or NPDR models are refit. The implementation was validated by exact reproduction of
the primary K = 50 results before application to K = 25 and K = 100. The main qualitative patterns were broadly consistent
across thresholds, although the magnitude of feature overlap and rank displacement varied with K.

Bootstrap stability was retained at K = 50 because the original bootstrap analysis stored replicate-level stability
summaries rather than complete feature rankings, and evaluation at additional K values would require rerunning the full
resampling procedure.

## Input files (not included)

The scripts expect the following files in the data directory (see below):

* `gene_with_label.csv` – real data: one row per subject, 1,783 gene columns plus `class` (HC/MDD), 157 subjects;
* `synthetic_rna_arf.csv`, `synthetic_ctgan.csv`, `synthetic_rna_tvae.csv` – synthetic data sets generated from the real data
  (200 subjects each).

## Running the scripts

All scripts read and write below one **data directory**, set with the environment variable `SYNTH_EVAL_BASE_DIR`
(default: the current directory); some folders accept their own variable (`HOLDOUT_BASE_DIR`, `RD_BASE_DIR`,
`PRIV_BASE_DIR`, `KS_DIR`), which takes precedence. Example:

```bash
export SYNTH_EVAL_BASE_DIR="$HOME/synth_eval_data"      # folder with the input files above (kept outside the repository)
Rscript revision_analyses/R6_privacy_fold_auc/09_privacy_ridge_domain_separator.R
```

Results are written below `$SYNTH_EVAL_BASE_DIR/reviewer_revision/...`; the `.gitignore` of this repository excludes these folders.

## Data availability

The real gene-expression dataset is not included in this repository because of data-use restrictions, and the synthetic
datasets generated from it are not shared either. The repository provides the preprocessing, synthetic-data generation,
evaluation and simulation code needed to reproduce the computational workflow for researchers with authorized access to the
underlying real data.

## Notes on reproducibility

* Unless otherwise specified, the analyses use a base random seed of 42, with deterministic seed offsets applied across
  repeated model fits and resampling procedures. Seeds for the simulated benchmark are in
  `revision_analyses/R5_known_signal/seeds.csv`.
* Synthetic-data generation is stochastic: regenerated synthetic datasets are valid but are not expected to be byte-identical
  to the datasets used for the manuscript. Because the synthetic data are not shared, the exact numerical results can be
  reproduced only with access to the real data and the archived synthetic files.
* CTGAN and TVAE fits take roughly 30–45 minutes (CTGAN, 900 epochs) on a laptop CPU; the repeated-draw analyses take several
  hours in total.

## Software

The primary evaluation pipeline uses R; CTGAN, TVAE and the quantum–adversarial extension use Python.

* R packages: see `requirements_R.txt`; Python packages: see `requirements_python.txt`.
* `Rscript session_info.R` writes `session_info_R.txt` (R and package versions). Run it in the installation used for the
  analyses and keep the file with the results.
* Model hyperparameters, seeds and evaluation parameters are reported in the Methods of the manuscript and in the scripts.


## Contact

Setareh Akbari (setareh-akbari@utulsa.edu).
  
