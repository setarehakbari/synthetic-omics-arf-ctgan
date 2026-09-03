############################################################
# ARF_generate.R
# Clean ARF synthetic-data generation script for the manuscript
#
# Purpose:
#   1) Load the canonical preprocessed real gene-expression data
#   2) Fit an Adversarial Random Forest (ARF)
#   3) Estimate the ARF density with FORDE
#   4) Generate 200 synthetic samples with FORGE
#   5) Validate and save the regenerated synthetic dataset
#
# IMPORTANT:
#   - This script NEVER overwrites the archived manuscript file
#     'synthetic_rna_arf.csv'.
#   - Regenerated data are written to
#     'synthetic_rna_arf_regenerated.csv'.
#   - Regenerated observations are not expected to be byte-identical
#     to the archived manuscript dataset unless the full historical
#     software environment and random-state behavior are reproduced.
############################################################

## ========== 0) User settings ==========
SEED <- 42
N_SYNTH <- 200
ARF_NUM_TREES <- 100
FORDE_ALPHA <- 0.1
FORDE_EPSILON <- 0.1

base_dir <- path.expand("~/Desktop/My paper")
input_csv <- file.path(base_dir, "gene_with_label.csv")

# Archived dataset used for the manuscript analyses.
# This file is READ-ONLY from the perspective of this script.
archived_csv <- file.path(base_dir, "synthetic_rna_arf.csv")

# Safe regenerated output: intentionally different filename.
output_csv <- file.path(base_dir, "synthetic_rna_arf_regenerated.csv")

# Auxiliary reproducibility information.
aux_dir <- file.path(base_dir, "ARF_generation")
dir.create(aux_dir, showWarnings = FALSE, recursive = TRUE)
summary_csv <- file.path(aux_dir, "arf_generation_summary.csv")
session_txt <- file.path(aux_dir, "arf_sessionInfo.txt")

# Safety switch: by default, do not overwrite a previous regenerated run.
ALLOW_OVERWRITE_REGENERATED <- FALSE

## ========== 1) Package check ==========
if (!requireNamespace("arf", quietly = TRUE)) {
  stop(
    "Package 'arf' is not installed. Install the ARF package used for the study, ",
    "then re-run this script."
  )
}

suppressPackageStartupMessages(library(arf))

## ========== 2) Safety checks ==========
if (!file.exists(input_csv)) {
  stop("Input file not found: ", input_csv)
}

if (normalizePath(output_csv, mustWork = FALSE) ==
    normalizePath(archived_csv, mustWork = FALSE)) {
  stop("Safety stop: regenerated output cannot be the archived manuscript file.")
}

if (file.exists(output_csv) && !ALLOW_OVERWRITE_REGENERATED) {
  stop(
    "Regenerated output already exists: ", output_csv, "\n",
    "Rename/move the existing regenerated file, or set ",
    "ALLOW_OVERWRITE_REGENERATED <- TRUE if you intentionally want to replace it."
  )
}

## ========== 3) Load canonical real data ==========
cat("Loading canonical real dataset ...\n")
real <- read.csv(
  input_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"class" %in% names(real)) {
  stop("Input data must contain a column named 'class'.")
}

# Manuscript dataset checks: 157 samples, 1783 gene features + class.
EXPECTED_REAL_ROWS <- 157
EXPECTED_GENE_FEATURES <- 1783

genes <- setdiff(names(real), "class")

if (nrow(real) != EXPECTED_REAL_ROWS) {
  stop(
    "Unexpected number of real samples. Expected ", EXPECTED_REAL_ROWS,
    ", found ", nrow(real), "."
  )
}

if (length(genes) != EXPECTED_GENE_FEATURES) {
  stop(
    "Unexpected number of gene features. Expected ", EXPECTED_GENE_FEATURES,
    ", found ", length(genes), "."
  )
}

# Validate class labels and keep the original 0/1 coding as a factor for ARF.
class_chr <- as.character(real$class)
if (!setequal(unique(class_chr), c("0", "1"))) {
  stop("The 'class' column must contain only 0 and 1.")
}
real$class <- factor(class_chr, levels = c("0", "1"))

# Gene features must be numeric and finite.
for (g in genes) {
  real[[g]] <- suppressWarnings(as.numeric(real[[g]]))
}

if (anyNA(real[, genes, drop = FALSE])) {
  stop("Missing or non-numeric values were detected in the gene-expression matrix.")
}

real_matrix <- as.matrix(real[, genes, drop = FALSE])
if (any(!is.finite(real_matrix))) {
  stop("Non-finite values (Inf/-Inf) were detected in the gene-expression matrix.")
}

zero_var <- genes[vapply(real[, genes, drop = FALSE], function(x) sd(x) == 0, logical(1))]
if (length(zero_var) > 0) {
  stop("Zero-variance genes detected: ", paste(zero_var, collapse = ", "))
}

cat("Real shape:", nrow(real), "x", ncol(real), "\n")
cat("Real class distribution:\n")
print(table(real$class))

## ========== 4) Fit ARF ==========
set.seed(SEED)
cat("\nFitting Adversarial Random Forest ...\n")
t0 <- proc.time()[["elapsed"]]

arf_model <- adversarial_rf(
  real,
  num_trees = ARF_NUM_TREES
)

## ========== 5) Estimate density with FORDE ==========
cat("Estimating density with FORDE ...\n")
psi <- forde(
  arf_model,
  real,
  alpha = FORDE_ALPHA,
  epsilon = FORDE_EPSILON
)

## ========== 6) Generate synthetic data with FORGE ==========
cat("Generating", N_SYNTH, "synthetic samples with FORGE ...\n")
synthetic <- forge(
  psi,
  n_synth = N_SYNTH
)

training_and_generation_seconds <- proc.time()[["elapsed"]] - t0

synthetic <- as.data.frame(synthetic, check.names = FALSE)

## ========== 7) Validate regenerated synthetic data ==========
if (nrow(synthetic) != N_SYNTH) {
  stop("Unexpected number of synthetic rows. Expected ", N_SYNTH,
       ", found ", nrow(synthetic), ".")
}

if (!setequal(names(synthetic), names(real))) {
  missing_cols <- setdiff(names(real), names(synthetic))
  extra_cols <- setdiff(names(synthetic), names(real))
  stop(
    "Synthetic columns do not match the real-data columns.\n",
    "Missing: ", paste(missing_cols, collapse = ", "), "\n",
    "Extra: ", paste(extra_cols, collapse = ", ")
  )
}

# Enforce the exact canonical column order.
synthetic <- synthetic[, names(real), drop = FALSE]

# Standardize class back to numeric 0/1 for downstream scripts.
synthetic$class <- suppressWarnings(as.integer(as.character(synthetic$class)))
if (anyNA(synthetic$class) || !all(synthetic$class %in% c(0L, 1L))) {
  stop("Synthetic 'class' values are not valid 0/1 labels.")
}

# Force gene columns to numeric and validate.
for (g in genes) {
  synthetic[[g]] <- suppressWarnings(as.numeric(synthetic[[g]]))
}

if (anyNA(synthetic[, genes, drop = FALSE])) {
  stop("Missing/non-numeric values detected in regenerated synthetic genes.")
}

synthetic_matrix <- as.matrix(synthetic[, genes, drop = FALSE])
if (any(!is.finite(synthetic_matrix))) {
  stop("Non-finite values detected in regenerated synthetic genes.")
}

if (ncol(synthetic) != EXPECTED_GENE_FEATURES + 1L) {
  stop("Unexpected synthetic column count: ", ncol(synthetic))
}

cat("Synthetic shape:", nrow(synthetic), "x", ncol(synthetic), "\n")
cat("Synthetic class distribution:\n")
print(table(synthetic$class))
cat("Missing values:", sum(is.na(synthetic)), "\n")
cat("Duplicated rows:", sum(duplicated(synthetic)), "\n")

## ========== 8) Save regenerated data ==========
write.csv(
  synthetic,
  output_csv,
  row.names = FALSE,
  quote = FALSE
)
cat("\nRegenerated synthetic data saved to:\n", output_csv, "\n", sep = "")

# Explicitly report the archived file status without modifying it.
if (file.exists(archived_csv)) {
  cat("Archived manuscript dataset preserved unchanged at:\n",
      archived_csv, "\n", sep = "")
} else {
  cat("Note: archived manuscript dataset was not found at:\n",
      archived_csv, "\n", sep = "")
}

## ========== 9) Save reproducibility summary ==========
real_counts <- table(real$class)
syn_counts <- table(factor(synthetic$class, levels = c(0, 1)))

summary_df <- data.frame(
  input_file = input_csv,
  archived_manuscript_file = archived_csv,
  regenerated_output_file = output_csv,
  random_seed = SEED,
  real_rows = nrow(real),
  gene_features = length(genes),
  real_class_0 = as.integer(real_counts["0"]),
  real_class_1 = as.integer(real_counts["1"]),
  synthetic_rows = nrow(synthetic),
  synthetic_class_0 = as.integer(syn_counts["0"]),
  synthetic_class_1 = as.integer(syn_counts["1"]),
  arf_num_trees = ARF_NUM_TREES,
  forde_alpha = FORDE_ALPHA,
  forde_epsilon = FORDE_EPSILON,
  training_and_generation_seconds = training_and_generation_seconds,
  duplicated_synthetic_rows = sum(duplicated(synthetic)),
  arf_package_version = as.character(utils::packageVersion("arf")),
  R_version = R.version.string,
  platform = R.version$platform,
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_csv, row.names = FALSE, quote = TRUE)
capture.output(sessionInfo(), file = session_txt)

cat("Generation summary saved to:\n", summary_csv, "\n", sep = "")
cat("R session information saved to:\n", session_txt, "\n", sep = "")

cat("\nARF generation completed successfully.\n")
############################################################
# End of ARF_generate.R
############################################################
