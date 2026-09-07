# preprocess_real_data.R
# Reproduce the real-data preprocessing used in the manuscript:
# quantile normalization -> log2 transformation -> CoV filtering -> zero-variance removal
# -> sample-by-gene matrix -> binary class label.
#
# NOTE:
# - Raw/real data are not distributed with this repository because of data-use restrictions.
# - Users with authorized access should update the two input paths below.
# - This script intentionally reproduces the preprocessing steps used for the manuscript dataset.

suppressPackageStartupMessages({
  library(dplyr)
  library(preprocessCore)
})

# -------------------------------------------------------------------------
# 1. User-configurable paths
# -------------------------------------------------------------------------

EXPRESSION_RDATA <- "path/to/sense.filtered.cpm.Rdata"
PHENOTYPE_CSV    <- "path/to/Demographic_symptom.csv"
OUTPUT_CSV       <- "gene_with_label.csv"

# Expected structure of the final manuscript dataset
EXPECTED_SAMPLES <- 157L
EXPECTED_GENES   <- 1783L

# Coefficient-of-variation threshold used in the manuscript
COV_THRESHOLD <- 0.025

# -------------------------------------------------------------------------
# 2. Load expression data
# -------------------------------------------------------------------------

load(EXPRESSION_RDATA)

if (!exists("sense.filtered.cpm")) {
  stop(
    "The RData file must contain an object named 'sense.filtered.cpm'."
  )
}

expr <- as.matrix(sense.filtered.cpm)

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene row names and sample column names.")
}

# -------------------------------------------------------------------------
# 3. Load and align phenotype data
# -------------------------------------------------------------------------

subject.attrs <- read.csv(
  PHENOTYPE_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_pheno_cols <- c("X", "Diag")
missing_pheno_cols <- setdiff(required_pheno_cols, names(subject.attrs))
if (length(missing_pheno_cols) > 0) {
  stop(
    "Phenotype file is missing required column(s): ",
    paste(missing_pheno_cols, collapse = ", ")
  )
}

# Keep phenotype records corresponding to expression samples.
phenos.df <- subject.attrs %>%
  filter(X %in% colnames(expr)) %>%
  dplyr::select(X, Diag)

if (anyDuplicated(phenos.df$X)) {
  stop("Duplicate subject IDs were found in the phenotype file.")
}

# Explicitly align phenotype rows to the expression-matrix column order.
match_idx <- match(colnames(expr), phenos.df$X)

if (anyNA(match_idx)) {
  missing_ids <- colnames(expr)[is.na(match_idx)]
  stop(
    "No phenotype record was found for expression sample(s): ",
    paste(missing_ids, collapse = ", ")
  )
}

phenos.df <- phenos.df[match_idx, , drop = FALSE]

if (!identical(as.character(phenos.df$X), colnames(expr))) {
  stop("Phenotype and expression sample ordering could not be aligned.")
}

# -------------------------------------------------------------------------
# 4. Quantile normalization and log2 transformation
# -------------------------------------------------------------------------

expr_quantile <- preprocessCore::normalize.quantiles(expr)
expr_log2 <- log2(expr_quantile)

colnames(expr_log2) <- colnames(expr)
rownames(expr_log2) <- rownames(expr)

if (any(!is.finite(expr_log2))) {
  stop(
    "Non-finite values were produced after log2 transformation. ",
    "Check the authorized source expression data."
  )
}

# -------------------------------------------------------------------------
# 5. Coefficient-of-variation and zero-variance filtering
# -------------------------------------------------------------------------

# CoV = sd(x) / abs(mean(x)), matching the manuscript preprocessing.
cov_values <- apply(
  expr_log2,
  1,
  function(x) sd(x) / abs(mean(x))
)

sd_values <- apply(
  expr_log2,
  1,
  sd
)

keep <- is.finite(cov_values) &
        cov_values < COV_THRESHOLD &
        sd_values > 0

expr_filtered <- expr_log2[keep, , drop = FALSE]

cat("Genes before filtering:", nrow(expr_log2), "\n")
cat("Genes retained after CoV and zero-variance filtering:",
    nrow(expr_filtered), "\n")

# -------------------------------------------------------------------------
# 6. Construct sample-by-gene analysis matrix and class label
# -------------------------------------------------------------------------

real <- as.data.frame(
  t(expr_filtered),
  check.names = FALSE
)

diag_chr <- as.character(phenos.df$Diag)

if (!all(diag_chr %in% c("HC", "MDD"))) {
  stop(
    "Unexpected diagnostic labels. Expected only 'HC' and 'MDD'; found: ",
    paste(sort(unique(diag_chr)), collapse = ", ")
  )
}

# HC = 0, MDD = 1, as used in the manuscript analysis.
real$class <- ifelse(diag_chr == "MDD", 1L, 0L)

# -------------------------------------------------------------------------
# 7. Validation
# -------------------------------------------------------------------------

if (nrow(real) != EXPECTED_SAMPLES) {
  stop(
    "Unexpected number of samples. Expected ", EXPECTED_SAMPLES,
    ", found ", nrow(real), "."
  )
}

gene_cols <- setdiff(names(real), "class")

if (length(gene_cols) != EXPECTED_GENES) {
  stop(
    "Unexpected number of retained genes. Expected ", EXPECTED_GENES,
    ", found ", length(gene_cols), "."
  )
}

if (anyNA(real)) {
  stop("Missing values detected in the final preprocessed dataset.")
}

gene_matrix <- as.matrix(real[, gene_cols, drop = FALSE])

if (any(!is.finite(gene_matrix))) {
  stop("Non-finite values detected in the final preprocessed gene matrix.")
}

zero_var_final <- gene_cols[
  vapply(real[, gene_cols, drop = FALSE], function(x) sd(x) == 0, logical(1))
]

if (length(zero_var_final) > 0) {
  stop(
    "Zero-variance genes remained after preprocessing: ",
    paste(zero_var_final, collapse = ", ")
  )
}

if (!setequal(unique(real$class), c(0L, 1L))) {
  stop("Final class column does not contain both 0 and 1 labels.")
}

cat("\nFinal dataset:", nrow(real), "samples x",
    length(gene_cols), "genes + class label\n")
cat("Class distribution:\n")
print(table(real$class))

# -------------------------------------------------------------------------
# 8. Save
# -------------------------------------------------------------------------

write.csv(
  real,
  OUTPUT_CSV,
  row.names = FALSE,
  quote = TRUE
)

cat("\nSaved preprocessed dataset to:", OUTPUT_CSV, "\n")
