############################################################
# 07_privacy_fold_auc.R
# Domain-separator AUC computed within each validation fold
# (5 x 5 repeated stratified CV) with a FIXED positive class
# (real = 1, synthetic = 0), for ARF, CTGAN and TVAE.
#
# Same design, seeds and settings as section 8 of
# evaluate_all_generators.R, but ONLY this analysis (a few minutes).
#
# Run (Terminal):  Rscript 07_privacy_fold_auc.R
# Needs in the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory) : gene_with_label.csv,
#   synthetic_rna_arf.csv, synthetic_ctgan.csv, synthetic_rna_tvae.csv
# Writes to the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory)/reviewer_revision/R6_privacy_fold_auc/ :
#   04_privacy_fold_auc.csv  (one row per fold and generator)
#   04_privacy_summary.csv   (mean, SD, min, max of fold AUC; pooled AUC)
############################################################
suppressPackageStartupMessages({ library(tidyverse); library(ranger); library(pROC); library(rsample) })
SEED <- 42; PRIVACY_TREES <- 2000; PRIVACY_FOLDS <- 5; PRIVACY_REPEATS <- 5
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

gen_names <- c("ARF", "CTGAN", "TVAE"); fold_all <- list(); summ_all <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; cat("Real vs", g, "...\n")
  dom <- bind_rows(mutate(real, source = "Real"), mutate(syn, source = "Synthetic"))
  dom$source <- factor(dom$source, levels = c("Synthetic", "Real"))   # Real = positive class
  set.seed(SEED + 7000 + g_idx)
  folds <- vfold_cv(dom, v = PRIVACY_FOLDS, repeats = PRIVACY_REPEATS, strata = source)
  rows <- vector("list", nrow(folds)); oof <- vector("list", nrow(folds))
  for (i in seq_len(nrow(folds))) {
    train <- analysis(folds$splits[[i]]); test <- assessment(folds$splits[[i]])
    m <- ranger(source ~ . - class, data = train, probability = TRUE, num.trees = PRIVACY_TREES,
                seed = SEED + 8000 + 1000 * g_idx + i)
    pr <- predict(m, data = test)$predictions[, "Real"]
    auc_i <- as.numeric(pROC::auc(pROC::roc(test$source, pr, levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)))
    rows[[i]] <- tibble(generator = g, fold_id = i, fold_auc = auc_i)
    oof[[i]] <- tibble(truth = test$source, prob_real = pr)
  }
  fa <- bind_rows(rows); o <- bind_rows(oof)
  pooled <- as.numeric(pROC::auc(pROC::roc(o$truth, o$prob_real, levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)))
  fold_all[[g]] <- fa
  summ_all[[g]] <- tibble(generator = g, n_folds = nrow(fa), mean_fold_auc = mean(fa$fold_auc), sd_fold_auc = sd(fa$fold_auc),
                          min_fold_auc = min(fa$fold_auc), max_fold_auc = max(fa$fold_auc),
                          pooled_auc = pooled, normalized_auc = max(pooled, 1 - pooled))
}
write_csv(bind_rows(fold_all), file.path(out_dir, "04_privacy_fold_auc.csv"))
write_csv(bind_rows(summ_all), file.path(out_dir, "04_privacy_summary.csv"))
print(bind_rows(summ_all), width = Inf)
cat("Results written to:", out_dir, "\n")
