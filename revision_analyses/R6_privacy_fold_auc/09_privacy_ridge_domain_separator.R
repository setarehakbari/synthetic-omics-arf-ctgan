############################################################
# 09_privacy_ridge_domain_separator.R
# Domain-separator AUC with a NON-TREE classifier (ridge logistic regression).
#
# Reviewer 4 (comment 2): the random-forest domain separator is tree-based, as is ARF.
# This script repeats 07_privacy_fold_auc.R with ridge logistic regression instead of a
# random forest: same data, same repeated stratified 5 x 5 cross-validation, SAME folds
# (same seeds), real = positive class (real = 1, synthetic = 0), AUC computed within each of
# the 25 validation folds for ARF, CTGAN and TVAE.
#
# Ridge classifier: predictors standardized with the training fold only; penalty chosen by
# stratified 5-fold cross-validation inside the training fold (lambda.min).
#
# Run (Terminal):  Rscript 09_privacy_ridge_domain_separator.R      (a few minutes)
# Needs in the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory) : gene_with_label.csv, synthetic_rna_arf.csv,
#   synthetic_ctgan.csv, synthetic_rna_tvae.csv
# Writes to the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory)/reviewer_revision/R6_privacy_fold_auc/ :
#   04b_privacy_ridge_fold_auc.csv  (one row per fold and generator)
#   04b_privacy_ridge_summary.csv   (mean, SD, min, max of fold AUC; pooled AUC)
############################################################
suppressPackageStartupMessages({ library(tidyverse); library(glmnet); library(pROC); library(rsample) })
SEED <- 42; PRIVACY_FOLDS <- 5; PRIVACY_REPEATS <- 5; INNER_FOLDS <- 5
base_dir <- Sys.getenv("PRIV_BASE_DIR", Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
out_dir  <- file.path(base_dir, "reviewer_revision", "R6_privacy_fold_auc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
paths <- c(Real = "gene_with_label.csv", ARF = "synthetic_rna_arf.csv",
           CTGAN = "synthetic_ctgan.csv", TVAE = "synthetic_rna_tvae.csv")
paths <- setNames(file.path(base_dir, paths), names(paths))
if (!all(file.exists(paths))) stop("Missing: ", paste(paths[!file.exists(paths)], collapse = ", "))

fix_class <- function(x) {
  z <- trimws(as.character(x))
  if (all(unique(z) %in% c("0", "1"))) return(factor(z, levels = c("0", "1"), labels = c("HC", "MDD")))
  z <- toupper(z); stopifnot(all(unique(z) %in% c("HC", "MDD"))); factor(z, levels = c("HC", "MDD"))
}
real <- read.csv(paths["Real"], check.names = FALSE, stringsAsFactors = FALSE)
genes <- setdiff(names(real), "class"); stopifnot(nrow(real) == 157, length(genes) == 1783)
real$class <- fix_class(real$class)
load_syn <- function(p) {
  s <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(setequal(setdiff(names(s), "class"), genes), nrow(s) == 200)
  s <- s[, c(genes, "class")]; s$class <- fix_class(s$class); s
}
datasets <- list(ARF = load_syn(paths["ARF"]), CTGAN = load_syn(paths["CTGAN"]), TVAE = load_syn(paths["TVAE"]))

ridge_scores <- function(train, test, seed) {
  xtr <- as.matrix(train[, genes]); mu <- colMeans(xtr); s <- apply(xtr, 2, sd); s[s == 0] <- 1
  z <- function(x) sweep(sweep(as.matrix(x[, genes]), 2, mu), 2, s, "/")
  y <- as.integer(train$source == "Real")
  set.seed(seed); fid <- integer(nrow(train))
  for (lv in c(0, 1)) { ix <- which(y == lv); fid[ix] <- sample(rep(seq_len(INNER_FOLDS), length.out = length(ix))) }
  cv <- cv.glmnet(z(train), y, family = "binomial", alpha = 0, standardize = FALSE, foldid = fid)
  out <- as.numeric(predict(cv, newx = z(test), s = "lambda.min", type = "response"))
  attr(out, "at_max") <- isTRUE(cv$lambda.min >= max(cv$lambda))   # TRUE = largest penalty on the glmnet path was selected
  attr(out, "n_unique") <- length(unique(out))                      # 1 = every test sample received the same score (AUC = 0.5)
  attr(out, "max_abs_coef") <- max(abs(as.numeric(coef(cv, s = "lambda.min"))[-1]))
  out
}

gen_names <- c("ARF", "CTGAN", "TVAE"); fold_all <- list(); summ_all <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; cat("Real vs", g, "(ridge) ...\n")
  dom <- bind_rows(mutate(real, source = "Real"), mutate(syn, source = "Synthetic"))
  dom$source <- factor(dom$source, levels = c("Synthetic", "Real"))   # Real = positive class
  set.seed(SEED + 7000 + g_idx)                                       # same folds as 07_privacy_fold_auc.R
  folds <- vfold_cv(dom, v = PRIVACY_FOLDS, repeats = PRIVACY_REPEATS, strata = source)
  rows <- vector("list", nrow(folds)); oof <- vector("list", nrow(folds))
  for (i in seq_len(nrow(folds))) {
    train <- analysis(folds$splits[[i]]); test <- assessment(folds$splits[[i]])
    pr <- ridge_scores(train, test, SEED + 9000 + 1000 * g_idx + i)
    auc_i <- as.numeric(pROC::auc(pROC::roc(test$source, pr, levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)))
    rows[[i]] <- tibble(generator = g, fold_id = i, fold_auc = auc_i, intercept_only = attr(pr, "at_max"), n_unique_pred = attr(pr, "n_unique"), max_abs_coef = attr(pr, "max_abs_coef"))
    oof[[i]] <- tibble(truth = test$source, prob_real = pr)
    if (i %% 5 == 0) cat(sprintf("  fold %d/%d\n", i, nrow(folds)))
  }
  fa <- bind_rows(rows); o <- bind_rows(oof)
  pooled <- as.numeric(pROC::auc(pROC::roc(o$truth, o$prob_real, levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)))
  fold_all[[g]] <- fa
  summ_all[[g]] <- tibble(generator = g, n_folds = nrow(fa), mean_fold_auc = mean(fa$fold_auc), sd_fold_auc = sd(fa$fold_auc),
                          min_fold_auc = min(fa$fold_auc), max_fold_auc = max(fa$fold_auc), n_folds_intercept_only = sum(fa$intercept_only), n_folds_constant_scores = sum(fa$n_unique_pred == 1),
                          pooled_auc = pooled, normalized_auc = max(pooled, 1 - pooled))
}
write_csv(bind_rows(fold_all), file.path(out_dir, "04b_privacy_ridge_fold_auc.csv"))
write_csv(bind_rows(summ_all), file.path(out_dir, "04b_privacy_ridge_summary.csv"))
print(bind_rows(summ_all), width = Inf)
cat("Results written to:", out_dir, "\n")
