############################################################
# 03_evaluate_holdout.R   (stage 3 of 3)
#
# Held-out evaluation of synthetic augmentation.
# Addresses Reviewer 3 Major Comment 2 and Reviewer 4 points 1-3.
#
# For every split (stage 1) and every generator that has a
# synthetic file (stage 2), using the SAME held-out real subjects:
#
#  A) Utility, scored on held-out REAL subjects only
#       training sets : Real-train | Synthetic only | Real-train + Synthetic
#       classifiers   : Random forest (tree-based)
#                       penalised logistic regression (ridge, alpha = 0 by default; non-tree;
#                       standardisation and penalty chosen by inner
#                       cross-validation inside the training set only)
#       metrics       : accuracy, balanced accuracy, AUC
#                       (positive class fixed to MDD; AUC is not
#                       folded with max(AUC, 1 - AUC))
#
#  B) Dilution of the real-train Top-K ranking in every split
#       Observed          : Real-train vs Real-train + Synthetic
#       Control A         : Real-train + rows resampled from Real-train
#                           (same class counts as the synthetic set)
#       Control B         : Real-train re-ranked with another RF seed
#     Every comparison is repeated N_REPS times inside each split
#     (Observed: different RF seeds on the same real + synthetic data;
#     Control A: different resampled rows and RF seeds; Control B:
#     different RF seeds). Repetition 1 uses the seed stored in
#     seeds.csv; repetitions 2..N_REPS use seeds derived from it.
#     The repetitions measure ranking noise given the training data;
#     they do not repeat the synthetic-data generation itself.
#
# Everything is descriptive: summaries across splits report mean,
# SD, min and max. Splits share subjects, so they are not
# independent replicates and the ranges are not confidence intervals.
#
# Usage
#   Rscript 03_evaluate_holdout.R
#   HOLDOUT_SPLITS=1,2 Rscript 03_evaluate_holdout.R   # selected splits
#   HOLDOUT_QUICK=1 Rscript 03_evaluate_holdout.R      # fewer trees: smoke test only
# Finished splits are kept unless HOLDOUT_OVERWRITE=1.
############################################################

## ========== 0) Settings ==========
base_dir <- Sys.getenv("HOLDOUT_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
out_dir  <- file.path(base_dir, "reviewer_revision", "R3_major2_holdout")
QUICK     <- identical(Sys.getenv("HOLDOUT_QUICK", "0"), "1")
OVERWRITE <- identical(Sys.getenv("HOLDOUT_OVERWRITE", "0"), "1")
SFX       <- if (QUICK) "_QUICK" else ""   # quick smoke tests never touch the file names of the real run
EVAL_TREES <- if (QUICK) 200L  else 2000L   # random-forest classifier
RANK_TREES <- if (QUICK) 300L  else 3000L   # permutation-importance ranking (as in the manuscript)
TOP_K      <- 50L
N_REPS     <- if (QUICK) 3L else as.integer(Sys.getenv("HOLDOUT_CTRL_REPS", "10"))  # repetitions of every ranking comparison
LR_ALPHA   <- as.numeric(Sys.getenv("HOLDOUT_LR_ALPHA", "0"))  # glmnet mixing: 0 = ridge (default), 0.5 = elastic net
INNER_FOLDS <- 5L
GENERATORS <- c("ARF", "CTGAN", "CTGAN_long", "TVAE")

suppressPackageStartupMessages({
  library(ranger); library(glmnet); library(pROC); library(dplyr); library(tibble)
})
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

seeds <- read.csv(file.path(out_dir, "seeds.csv"), stringsAsFactors = FALSE)
get_seed <- function(i, role) {
  s <- seeds$seed[seeds$split == i & seeds$role == role]
  if (length(s) != 1) stop("Seed not found for split ", i, ", role ", role)
  s
}

# Repetition 1 = the stored seed; repetitions 2..n are drawn (without replacement) from a generator started by it.
derive_seeds <- function(base, n) {
  if (n <= 1L) return(base)
  set.seed(base)
  c(base, sample.int(2147483647L, n - 1L))
}

## ========== 1) Helpers ==========
as_class <-function(x) factor(ifelse(x == 1, "MDD", "HC"), levels = c("HC", "MDD"))

rf_fit_predict <- function(train, test, seed) {
  m <- ranger(class ~ ., data = train, num.trees = EVAL_TREES, probability = TRUE, seed = seed)
  predict(m, data = test)$predictions[, "MDD"]
}
lr_fit_predict <- function(train, test, seed) {
  genes <- setdiff(names(train), "class")
  set.seed(seed)
  cv <- cv.glmnet(as.matrix(train[, genes]), train$class, family = "binomial",
                  alpha = LR_ALPHA, nfolds = INNER_FOLDS, standardize = TRUE)
  as.numeric(predict(cv, newx = as.matrix(test[, genes]), s = "lambda.min", type = "response"))
}
score <- function(p, truth) {
  pred <- ifelse(p >= 0.5, "MDD", "HC")
  sens <- mean(pred[truth == "MDD"] == "MDD"); spec <- mean(pred[truth == "HC"] == "HC")
  auc  <- as.numeric(pROC::auc(pROC::roc(truth, p, levels = c("HC", "MDD"), direction = "<", quiet = TRUE)))
  c(accuracy = mean(pred == truth), balanced_accuracy = (sens + spec) / 2, auc = auc)
}
rf_rank <- function(df, seed) {
  m <- ranger(class ~ ., data = df, num.trees = RANK_TREES, importance = "permutation",
              scale.permutation.importance = TRUE, seed = seed)
  tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |> mutate(rank = row_number())
}
dilution_vs_real <- function(rank_real, rank_aug, K = TOP_K) {
  top <- rank_real |> filter(rank <= K) |> select(gene, rank_real = rank)
  d <- top |> left_join(select(rank_aug, gene, rank_aug = rank), by = "gene") |>
    mutate(abs_shift = abs(rank_aug - rank_real), dropped = rank_aug > K)
  top_aug <- rank_aug$gene[rank_aug$rank <= K]
  tibble(mean_abs_shift = mean(d$abs_shift), median_abs_shift = median(d$abs_shift),
         dropout_pct = 100 * mean(d$dropped),
         jaccard_top = length(intersect(top$gene, top_aug)) / length(union(top$gene, top_aug)))
}

## ========== 2) Per-split evaluation ==========
all_splits <- sort(as.integer(sub("split_", "", list.files(file.path(out_dir, "splits"), pattern = "^split_"))))
sel <- Sys.getenv("HOLDOUT_SPLITS", unset = "")
splits <- if (nzchar(sel)) as.integer(strsplit(sel, ",")[[1]]) else all_splits

for (i in splits) {
  d <- file.path(out_dir, "splits", sprintf("split_%02d", i))
  f_util <- file.path(d, sprintf("results_utility%s.csv", SFX)); f_dil <- file.path(d, sprintf("results_dilution%s.csv", SFX))
  # A split counts as finished only if (a) its results cover exactly the generators that have a synthetic
  # file now, and (b) no input (train/test/synthetic) is newer than the results.
  present <- GENERATORS[file.exists(file.path(d, sprintf("synthetic_%s.csv", GENERATORS)))]
  finished <- file.exists(f_util) && file.exists(f_dil) && {
    dd <- read.csv(f_dil); uu <- read.csv(f_util)
    done_g <- setdiff(unique(dd$generator), "none")
    inputs <- file.path(d, c("train.csv", "test.csv", sprintf("synthetic_%s.csv", present)))
    ("rep" %in% names(dd)) && max(dd$rep) == N_REPS && setequal(done_g, present) &&
      ("lr_alpha" %in% names(uu)) && all(uu$lr_alpha == LR_ALPHA) &&   # re-evaluate if the logistic-regression penalty changed
      all(file.mtime(inputs) < min(file.mtime(c(f_util, f_dil))))
  }
  if (finished && !OVERWRITE) { msg("split %02d: results up to date, skipping", i); next }
  if (file.exists(f_util) && !OVERWRITE) msg("split %02d: earlier results are incomplete or older than the inputs; re-evaluating", i)
  t_split <- proc.time()[["elapsed"]]

  train <- read.csv(file.path(d, "train.csv"), check.names = FALSE)
  test  <- read.csv(file.path(d, "test.csv"),  check.names = FALSE)
  genes <- setdiff(names(train), "class")
  stopifnot(identical(names(train), names(test)))
  train$class <- as_class(train$class); test$class <- as_class(test$class)
  truth <- test$class
  msg("\n=== split %02d: train %d | test %d (HC %d / MDD %d) | %d genes ===",
      i, nrow(train), nrow(test), sum(truth == "HC"), sum(truth == "MDD"), length(genes))

  util <- list(); dil <- list()
  add_util <- function(model, training_set, generator, p, n_rows) {
    util[[length(util) + 1]] <<- tibble(split = i, model = model, training_set = training_set, generator = generator,
                                        n_train_rows = n_rows, n_test = length(truth), lr_alpha = LR_ALPHA, !!!as.list(score(p, truth)))
  }
  evaluate_training_set <- function(tr, training_set, generator, cond_key) {
    for (model in c("RF", "LR")) {
      p <- tryCatch(
        if (model == "RF") rf_fit_predict(tr, test, get_seed(i, paste0("rf_", cond_key)))
        else               lr_fit_predict(tr, test, get_seed(i, paste0("lr_", cond_key))),
        error = function(e) { msg("   ! %s on %s failed: %s", model, cond_key, conditionMessage(e)); NULL })
      if (!is.null(p)) add_util(model, training_set, generator, p, nrow(tr))
    }
  }

  # --- real-only baseline and its reference ranking
  evaluate_training_set(train, "Real-train", "none", "Real")
  rank_real <- rf_rank(train, get_seed(i, "rank_Real"))

  # Seeds for the repeated rankings must not collide with each other or with any stored seed.
  derived_roles <- c("rank_ctrlB", paste0("rank_RealPlus_", GENERATORS), paste0("ctrlA_", GENERATORS),
                     paste0("rank_ctrlA_", GENERATORS))
  extra <- unlist(lapply(derived_roles, function(rl) derive_seeds(get_seed(i, rl), N_REPS)[-1]))
  if (anyDuplicated(c(seeds$seed, extra))) stop("Derived seed collision in split ", i, ".")
  add_dil <- function(comparison, generator, rep, rank_other) {
    dil[[length(dil) + 1]] <<- bind_cols(tibble(split = i, comparison = comparison, generator = generator, rep = rep),
                                         dilution_vs_real(rank_real, rank_other))
  }

  # --- Control B: same real-train data, different random-forest seeds
  sB <- derive_seeds(get_seed(i, "rank_ctrlB"), N_REPS)
  for (r in seq_len(N_REPS)) add_dil("Control B (RF seed only)", "none", r, rf_rank(train, sB[r]))

  # --- generators
  for (g in GENERATORS) {
    f_syn <- file.path(d, sprintf("synthetic_%s.csv", g))
    if (!file.exists(f_syn)) { msg("   %s: no synthetic file, skipped", g); next }
    syn <- read.csv(f_syn, check.names = FALSE)
    stopifnot(setequal(names(syn), names(train)))
    syn <- syn[, names(train), drop = FALSE]; syn$class <- as_class(syn$class)
    both <- rbind(train, syn)
    msg("   %s: synthetic %d rows (HC %d / MDD %d)", g, nrow(syn), sum(syn$class == "HC"), sum(syn$class == "MDD"))

    evaluate_training_set(syn,  "Synthetic only",         g, paste0("Syn_", g))
    evaluate_training_set(both, "Real-train + Synthetic", g, paste0("RealPlus_", g))

    # Observed dilution (same real + synthetic data, different RF seeds) and Control A
    # (real-train + rows resampled from real-train with the same class counts as the synthetic set;
    # a new resample and a new RF seed in every repetition)
    sObs   <- derive_seeds(get_seed(i, paste0("rank_RealPlus_", g)), N_REPS)
    sA_row <- derive_seeds(get_seed(i, paste0("ctrlA_", g)),         N_REPS)
    sA_rf  <- derive_seeds(get_seed(i, paste0("rank_ctrlA_", g)),    N_REPS)
    for (r in seq_len(N_REPS)) {
      add_dil("Observed (Real-train + Synthetic)", g, r, rf_rank(both, sObs[r]))
      set.seed(sA_row[r])
      idx <- c(sample(which(train$class == "HC"),  sum(syn$class == "HC"),  replace = TRUE),
               sample(which(train$class == "MDD"), sum(syn$class == "MDD"), replace = TRUE))
      add_dil("Control A (resampled Real-train)", g, r, rf_rank(rbind(train, train[idx, ]), sA_rf[r]))
    }
  }
  write.csv(bind_rows(util), f_util, row.names = FALSE)
  write.csv(bind_rows(dil),  f_dil,  row.names = FALSE)
  msg("split %02d finished in %.1f min", i, (proc.time()[["elapsed"]] - t_split) / 60)
}

## ========== 3) Summaries across the splits found on disk ==========
done <- all_splits[file.exists(file.path(out_dir, "splits", sprintf("split_%02d", all_splits), sprintf("results_utility%s.csv", SFX)))]
if (length(done) == 0) stop("No finished splits to summarise.")
read_all <- function(f) bind_rows(lapply(done, function(i) read.csv(file.path(out_dir, "splits", sprintf("split_%02d", i), f))))
util <- read_all(sprintf("results_utility%s.csv", SFX)); dil <- read_all(sprintf("results_dilution%s.csv", SFX))
write.csv(util, file.path(out_dir, sprintf("holdout_utility_per_split%s.csv", SFX)),  row.names = FALSE)
write.csv(dil,  file.path(out_dir, sprintf("holdout_dilution_per_split%s.csv", SFX)), row.names = FALSE)

# Completeness check: every generator should have been evaluated on the same splits
obs_rows <- dil[grepl("^Observed", dil$comparison), ]
ns <- vapply(GENERATORS, function(g) length(unique(obs_rows$split[obs_rows$generator == g])), integer(1))  # number of splits per generator
if (length(unique(as.integer(ns))) > 1 || any(ns == 0))
  msg("\nWARNING: generators were evaluated on different numbers of splits (%s). Finish all generators and re-run stage 3 before reporting results.",
      paste(names(ns), as.integer(ns), sep = "=", collapse = ", "))

metrics <- c("accuracy", "balanced_accuracy", "auc")
utility_summary <- util |>
  group_by(model, training_set, generator) |>
  summarise(n_splits = n(),
            across(all_of(metrics), list(mean = mean, sd = sd, min = min, max = max), .names = "{.col}_{.fn}"),
            .groups = "drop")
real_ref <- util |> filter(training_set == "Real-train") |>
  select(split, model, real_accuracy = accuracy, real_balanced_accuracy = balanced_accuracy, real_auc = auc)
delta <- util |> filter(training_set == "Real-train + Synthetic") |> left_join(real_ref, by = c("split", "model")) |>
  mutate(d_accuracy = accuracy - real_accuracy, d_balanced_accuracy = balanced_accuracy - real_balanced_accuracy,
         d_auc = auc - real_auc) |>
  group_by(model, generator) |>
  summarise(n_splits = n(), across(c(d_accuracy, d_balanced_accuracy, d_auc), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            n_splits_auc_higher_than_real = sum(d_auc > 0), .groups = "drop")
write.csv(utility_summary, file.path(out_dir, sprintf("holdout_utility_summary%s.csv", SFX)), row.names = FALSE)
write.csv(delta,           file.path(out_dir, sprintf("holdout_utility_delta_vs_real%s.csv", SFX)), row.names = FALSE)

dil_metrics <- c("mean_abs_shift", "median_abs_shift", "dropout_pct", "jaccard_top")
## Step 1: inside each split, summarise the repetitions (mean, min, max of the mean absolute rank shift).
rep_range <- dil |>
  group_by(split, comparison, generator) |>
  summarise(n_reps = n(), mean_shift = mean(mean_abs_shift), min_shift = min(mean_abs_shift),
            max_shift = max(mean_abs_shift), .groups = "drop")
## Step 2: across splits, summarise the split-level means (descriptive; splits share subjects).
dil_split <- dil |>
  group_by(split, comparison, generator) |>
  summarise(n_reps = n(), across(all_of(dil_metrics), mean), .groups = "drop")
dilution_summary <- dil_split |>
  group_by(comparison, generator) |>
  summarise(n_splits = n(), n_reps_per_split = max(n_reps),
            across(all_of(dil_metrics), list(mean = mean, sd = sd, min = min, max = max), .names = "{.col}_{.fn}"),
            .groups = "drop")
## Step 3: within each split, is the observed shift above the controls?
obs   <- rep_range |> filter(grepl("^Observed",  comparison)) |> select(split, generator, obs_mean = mean_shift, obs_min = min_shift)
ctrlA <- rep_range |> filter(grepl("^Control A", comparison)) |> select(split, generator, A_mean = mean_shift, A_max = max_shift)
ctrlB <- rep_range |> filter(grepl("^Control B", comparison)) |> select(split, B_mean = mean_shift, B_max = max_shift)
within_split <- obs |> left_join(ctrlA, by = c("split", "generator")) |> left_join(ctrlB, by = "split") |>
  group_by(generator) |>
  summarise(n_splits = n(),
            n_obs_mean_above_ctrlA_mean = sum(obs_mean > A_mean), n_obs_mean_above_ctrlB_mean = sum(obs_mean > B_mean),
            n_obs_min_above_ctrlA_max   = sum(obs_min  > A_max),  n_obs_min_above_ctrlB_max   = sum(obs_min  > B_max),
            .groups = "drop")
write.csv(dilution_summary, file.path(out_dir, sprintf("holdout_dilution_summary%s.csv", SFX)), row.names = FALSE)
write.csv(rep_range,        file.path(out_dir, sprintf("holdout_dilution_within_split_spread%s.csv", SFX)), row.names = FALSE)
write.csv(within_split,     file.path(out_dir, sprintf("holdout_dilution_within_split_comparison%s.csv", SFX)), row.names = FALSE)

# budget log written by stages 2a/2b
logs <- list.files(file.path(out_dir, "splits"), pattern = "^log_.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(logs) > 0) {
  lg <- bind_rows(lapply(logs, read.csv, stringsAsFactors = FALSE))
  write.csv(lg, file.path(out_dir, sprintf("generator_budget_log%s.csv", SFX)), row.names = FALSE)
}
writeLines(c(sprintf("splits summarised: %s", paste(done, collapse = ", ")),
             sprintf("EVAL_TREES=%d RANK_TREES=%d TOP_K=%d LR_ALPHA=%.2f INNER_FOLDS=%d QUICK=%s",
                     EVAL_TREES, RANK_TREES, TOP_K, LR_ALPHA, INNER_FOLDS, QUICK), "", capture.output(sessionInfo())),
           file.path(out_dir, sprintf("run_info%s.txt", SFX)))

msg("\n=============== UTILITY on held-out real subjects (mean over %d splits) ===============", length(done))
print(as.data.frame(utility_summary |> transmute(model, training_set, generator, n_splits,
      accuracy = round(accuracy_mean, 3), bal_acc = round(balanced_accuracy_mean, 3), auc = round(auc_mean, 3),
      auc_sd = round(auc_sd, 3))), row.names = FALSE)
options(width = 200)
msg("\n=============== DILUTION (mean over splits; each split = mean of %d repetitions) ===============", N_REPS)
print(as.data.frame(dilution_summary |> transmute(comparison, generator, n_splits, n_reps = n_reps_per_split,
      mean_abs_shift = round(mean_abs_shift_mean, 1), sd = round(mean_abs_shift_sd, 1),
      dropout_pct = round(dropout_pct_mean, 1), jaccard = round(jaccard_top_mean, 3))), row.names = FALSE)
msg("\n=============== Observed vs controls, split by split ===============")
print(as.data.frame(within_split), row.names = FALSE)
msg("\nStage 3 finished. Output: %s", out_dir)
