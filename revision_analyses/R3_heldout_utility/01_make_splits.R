############################################################
# 01_make_splits.R   (stage 1 of 3)
#
# Held-out real-data evaluation of synthetic augmentation.
# Addresses Reviewer 3 Major Comment 2 and Reviewer 4 points 2-3.
#
# What this script does
#   Builds repeated, stratified train/test splits of the real
#   subjects. For EVERY split, all preprocessing parameters are
#   estimated on the TRAINING subjects only and then applied,
#   unchanged, to the held-out TEST subjects:
#     (1) quantile-normalisation reference distribution
#     (2) log2 transform (genes with non-finite log2 in the
#         training data are dropped; no pseudo-count)
#     (3) coefficient-of-variation filter (0 < CoV < COV_MAX)
#         and removal of zero-variance genes
#   The synthetic-data generators (stage 2) then see TRAINING
#   subjects only, and all models (stage 3) are scored on the same
#   held-out real subjects.
#
# Preprocessing reproduced from the raw matrix
#   Running this script with VERIFY_FULL_DATA <- TRUE applies the
#   same steps to all 157 subjects and checks that the result is
#   numerically identical to the manuscript matrix
#   (gene_with_label.csv; 1,783 genes).
#
# Inputs  (set base_dir below, or the environment variable
#          HOLDOUT_BASE_DIR)
#   sense.filtered.cpm.Rdata   object `sense.filtered.cpm`
#                              (genes x subjects, CPM)
#   Demographic_symptom.csv    first column = subject ID,
#                              column `Diag` = HC / MDD
#   gene_with_label.csv        (optional, only for verification)
#
# Outputs (under <base_dir>/reviewer_revision/R3_major2_holdout/)
#   seeds.csv                  one collision-free seed per
#                              (split, role); used by stages 2-3
#   splits/split_XX/train.csv  samples x genes + class (0=HC,1=MDD)
#   splits/split_XX/test.csv   same genes, held-out subjects
#   split_membership.csv       which subject is train/test in
#                              which split
#   split_summary.csv          sizes, class counts, genes kept
#
# Data privacy: input and output files contain patient-level
# expression data and must NOT be committed to a public repository.
############################################################

## ========== 0) Settings ==========
base_dir <- Sys.getenv("HOLDOUT_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
raw_rdata     <- file.path(base_dir, "sense.filtered.cpm.Rdata")
demo_csv      <- file.path(base_dir, "Demographic_symptom.csv")
processed_csv <- file.path(base_dir, "gene_with_label.csv")  # optional (verification only)

## Check the inputs BEFORE anything is created on disk.
if (!file.exists(raw_rdata) || !file.exists(demo_csv))
  stop("Input files not found in: ", base_dir, "\n",
       "This folder must contain sense.filtered.cpm.Rdata and Demographic_symptom.csv.\n",
       "Set the folder first: in R  Sys.setenv(HOLDOUT_BASE_DIR = dirname(file.choose()))  ",
       "(then pick sense.filtered.cpm.Rdata), or in Terminal  export HOLDOUT_BASE_DIR=\"/path/to/folder\"")

SEED       <- 42L
N_SPLITS   <- as.integer(Sys.getenv("HOLDOUT_N_SPLITS", unset = "5"))
TEST_FRAC  <- 0.25     # fraction of each class held out (stratified)
COV_MAX    <- 0.025    # keep genes with 0 < CoV < COV_MAX (as in the manuscript)
VERIFY_FULL_DATA <- TRUE

out_dir    <- file.path(base_dir, "reviewer_revision", "R3_major2_holdout")
splits_dir <- file.path(out_dir, "splits")

## Safety guard: re-running this stage would overwrite splits and seeds while old synthetic data
## and results (made from the OLD splits) could stay on disk and be mixed with the new splits.
## Stage 1 therefore refuses to run if splits already exist. HOLDOUT_RESET=1 deletes ALL earlier
## stage 1-3 outputs in this folder first.
if (dir.exists(splits_dir) && length(list.files(splits_dir)) > 0) {
  if (!identical(Sys.getenv("HOLDOUT_RESET", "0"), "1"))
    stop("Splits already exist in ", splits_dir, ".\n",
         "Do not re-run stage 1 unless you want to start over. To delete ALL earlier outputs ",
         "(splits, synthetic data, results) and rebuild, run again with HOLDOUT_RESET=1.")
  unlink(out_dir, recursive = TRUE)
  cat("HOLDOUT_RESET=1: earlier outputs deleted.\n")
}
dir.create(splits_dir, recursive = TRUE, showWarnings = FALSE)

## Seed roles used by stages 2 and 3 (defined here, in one place).
GENERATORS <- c("ARF", "CTGAN", "CTGAN_long", "TVAE")
CONDITIONS <- c("Real", paste0("Syn_", GENERATORS), paste0("RealPlus_", GENERATORS))
SEED_ROLES <- c("split",
                paste0("gen_", GENERATORS),
                paste0("rf_",  CONDITIONS),
                paste0("lr_",  CONDITIONS),
                paste0("rank_", c("Real", paste0("RealPlus_", GENERATORS))),
                paste0("ctrlA_",      GENERATORS),   # resampling control: which rows are drawn
                paste0("rank_ctrlA_", GENERATORS),   # resampling control: RF ranking seed
                "rank_ctrlB")                        # seed-only control: RF ranking seed

msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ========== 1) Load raw data and labels ==========
stopifnot(file.exists(raw_rdata), file.exists(demo_csv))
e <- new.env()
load(raw_rdata, envir = e)
raw <- e[["sense.filtered.cpm"]]
stopifnot(is.matrix(raw), !anyNA(raw), !is.null(colnames(raw)), !is.null(rownames(raw)))

demo <- read.csv(demo_csv, check.names = FALSE, stringsAsFactors = FALSE)
ids  <- demo[[1]]
diag <- demo$Diag[match(colnames(raw), ids)]
if (anyNA(diag) || !all(diag %in% c("HC", "MDD")))
  stop("Every subject in the expression matrix needs Diag = HC or MDD in the demographic file.")
class01 <- ifelse(diag == "MDD", 1L, 0L)
names(class01) <- colnames(raw)
msg("Raw matrix: %d genes x %d subjects | HC = %d, MDD = %d",
    nrow(raw), ncol(raw), sum(class01 == 0), sum(class01 == 1))

## ========== 2) Preprocessing functions ==========
## (1) quantile normalisation: reference = row means of the column-sorted TRAIN matrix
qn_reference <- function(M) rowMeans(apply(M, 2, sort))
qn_apply <- function(M, ref) {
  out <- M
  for (j in seq_len(ncol(M))) {
    r <- rank(M[, j], ties.method = "average")
    out[, j] <- (ref[floor(r)] + ref[ceiling(r)]) / 2
  }
  out
}

## Fit all preprocessing parameters on `raw_train` (genes x train subjects).
preprocess_fit <- function(raw_train) {
  ref <- qn_reference(raw_train)
  l   <- suppressWarnings(log2(qn_apply(raw_train, ref)))
  finite <- apply(l, 1, function(v) all(is.finite(v)))
  l   <- l[finite, , drop = FALSE]
  s   <- apply(l, 1, sd)
  m   <- rowMeans(l)
  cv  <- s / abs(m)            # CoV on the log2 scale; all gene means are positive here, so abs() changes nothing
  keep <- is.finite(cv) & s > 0 & cv > 0 & cv < COV_MAX
  list(reference = ref, genes = rownames(l)[keep])
}

## Apply fitted parameters to any subjects (genes x subjects) -> subjects x genes.
preprocess_apply <- function(raw_any, fit) {
  l <- log2(qn_apply(raw_any, fit$reference)[fit$genes, , drop = FALSE])
  if (any(!is.finite(l))) stop("Non-finite values after preprocessing; check the reference distribution.")
  t(l)
}

## ========== 3) Optional check: reproduce the manuscript matrix ==========
if (VERIFY_FULL_DATA && file.exists(processed_csv)) {
  msg("\nVerification: preprocessing ALL subjects should reproduce the manuscript matrix ...")
  ms  <- read.csv(processed_csv, check.names = FALSE)
  ms_genes <- setdiff(names(ms), "class")
  fit_all <- preprocess_fit(raw)
  X_all   <- preprocess_apply(raw, fit_all)
  colnames(X_all) <- make.names(colnames(X_all), unique = TRUE)
  same_set <- setequal(colnames(X_all), ms_genes)
  max_diff <- if (same_set) max(abs(X_all[, ms_genes] - as.matrix(ms[, ms_genes]))) else NA_real_
  msg("  genes kept: %d (manuscript: %d) | identical gene set: %s | max |difference|: %s",
      ncol(X_all), length(ms_genes), same_set, format(max_diff, digits = 3))
  if (!same_set || !is.finite(max_diff) || max_diff > 1e-8)
    stop("Preprocessing does NOT reproduce the manuscript matrix; do not continue until this is resolved.")
  msg("  OK: manuscript preprocessing reproduced.")
} else if (VERIFY_FULL_DATA) {
  msg("Note: %s not found; skipping the reproduction check.", processed_csv)
}

## ========== 4) Collision-free seed bank ==========
set.seed(SEED + 777000L)
n_seeds <- N_SPLITS * length(SEED_ROLES)
bank <- sample.int(2147483647L, size = n_seeds)          # no replacement => no duplicates
stopifnot(!anyDuplicated(bank))
seeds <- data.frame(
  split = rep(seq_len(N_SPLITS), each = length(SEED_ROLES)),
  role  = rep(SEED_ROLES, times = N_SPLITS),
  seed  = bank,
  stringsAsFactors = FALSE
)
write.csv(seeds, file.path(out_dir, "seeds.csv"), row.names = FALSE)
get_seed <- function(i, role) seeds$seed[seeds$split == i & seeds$role == role]

## ========== 5) Stratified repeated splits ==========
membership <- list(); summary_rows <- list()
for (i in seq_len(N_SPLITS)) {
  set.seed(get_seed(i, "split"))
  test_ids <- unlist(lapply(c(0L, 1L), function(k) {
    pool <- names(class01)[class01 == k]
    sample(pool, size = round(TEST_FRAC * length(pool)))
  }))
  train_ids <- setdiff(colnames(raw), test_ids)

  fit   <- preprocess_fit(raw[, train_ids, drop = FALSE])
  Xtr   <- preprocess_apply(raw[, train_ids, drop = FALSE], fit)
  Xte   <- preprocess_apply(raw[, test_ids,  drop = FALSE], fit)
  colnames(Xtr) <- colnames(Xte) <- make.names(colnames(Xtr), unique = TRUE)

  tr <- data.frame(Xtr, class = class01[train_ids], check.names = FALSE, row.names = NULL)
  te <- data.frame(Xte, class = class01[test_ids],  check.names = FALSE, row.names = NULL)

  d <- file.path(splits_dir, sprintf("split_%02d", i))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  write.csv(tr, file.path(d, "train.csv"), row.names = FALSE)
  write.csv(te, file.path(d, "test.csv"),  row.names = FALSE)

  membership[[i]] <- data.frame(
    split = i, sample_id = c(train_ids, test_ids),
    set = rep(c("train", "test"), c(length(train_ids), length(test_ids))),
    class = c(class01[train_ids], class01[test_ids]))
  summary_rows[[i]] <- data.frame(
    split = i, n_train = length(train_ids), n_test = length(test_ids),
    train_HC = sum(tr$class == 0), train_MDD = sum(tr$class == 1),
    test_HC  = sum(te$class == 0), test_MDD  = sum(te$class == 1),
    genes_kept = ncol(Xtr))
  msg("split %02d: train %d (HC %d / MDD %d) | test %d (HC %d / MDD %d) | genes kept %d",
      i, nrow(tr), sum(tr$class == 0), sum(tr$class == 1),
      nrow(te), sum(te$class == 0), sum(te$class == 1), ncol(Xtr))
}
write.csv(do.call(rbind, membership),   file.path(out_dir, "split_membership.csv"), row.names = FALSE)
write.csv(do.call(rbind, summary_rows), file.path(out_dir, "split_summary.csv"),    row.names = FALSE)
msg("\nStage 1 finished. Output: %s", out_dir)
