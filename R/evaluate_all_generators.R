############################################################
# evaluate_all_generators.R
# Unified evaluation of ARF, CTGAN, and TVAE synthetic
# gene-expression datasets for the revised manuscript.
############################################################

## =========================================================
## 0) USER SETTINGS
## =========================================================
QUICK_TEST <- FALSE
SEED <- 42
TOP_K <- 50
TOP_BAR <- 20
UTILITY_REPS_FULL <- 100
UTILITY_TREES_FULL <- 5000
RANK_TREES_FULL <- 3000
PRIVACY_TREES_FULL <- 2000
PRIVACY_FOLDS <- 5
PRIVACY_REPEATS <- 5
STABILITY_B_FULL <- 50
STABILITY_TREES_FULL <- 2000
NPDR_METHOD <- "relieff"
NPDR_METRIC <- "manhattan"
NPDR_KNN <- 10
NPDR_PADJ <- "fdr"
UMAP_NEIGHBORS <- 15
UMAP_MIN_DIST <- 0.1

if (QUICK_TEST) {
  UTILITY_REPS <- 3
  UTILITY_TREES <- 300
  RANK_TREES <- 500
  PRIVACY_TREES <- 300
  STABILITY_B <- 3
  STABILITY_TREES <- 300
} else {
  UTILITY_REPS <- UTILITY_REPS_FULL
  UTILITY_TREES <- UTILITY_TREES_FULL
  RANK_TREES <- RANK_TREES_FULL
  PRIVACY_TREES <- PRIVACY_TREES_FULL
  STABILITY_B <- STABILITY_B_FULL
  STABILITY_TREES <- STABILITY_TREES_FULL
}

## =========================================================
## 1) PATHS
## =========================================================
base_dir <- path.expand("~/Desktop/My paper")
real_path <- file.path(base_dir, "gene_with_label.csv")
arf_path  <- file.path(base_dir, "synthetic_rna_arf.csv")
ctg_path  <- file.path(base_dir, "synthetic_ctgan.csv")
tvae_path <- file.path(base_dir, "synthetic_rna_tvae.csv")

out_dir <- file.path(base_dir, "reviewer_revision", "all_generators")
tables_dir  <- file.path(out_dir, "tables")
details_dir <- file.path(out_dir, "details")
figures_dir <- file.path(out_dir, "figures")
logs_dir    <- file.path(out_dir, "logs")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(details_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,    recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2) PACKAGES
## =========================================================
required_cran <- c("tidyverse", "ranger", "umap", "pROC", "rsample")
missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  stop(
    "Missing required CRAN package(s): ", paste(missing_cran, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_cran), collapse = ", "), "))"
  )
}
if (!requireNamespace("npdr", quietly = TRUE)) {
  stop(
    "Package 'npdr' is required.\nInstall it with:\n",
    "install.packages('remotes')\n",
    "remotes::install_github('insilico/npdr')"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(umap)
  library(pROC)
  library(rsample)
  library(npdr)
})
set.seed(SEED)

## =========================================================
## 3) HELPERS
## =========================================================
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  m <- mean(x); se <- sd(x) / sqrt(length(x))
  c(m - 1.96 * se, m + 1.96 * se)
}
jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}
fix_class <- function(x, dataset_name) {
  z <- trimws(as.character(x))
  uz <- unique(z)
  if (all(uz %in% c("0", "1"))) {
    return(factor(z, levels = c("0", "1"), labels = c("HC", "MDD")))
  }
  z_up <- toupper(z)
  if (all(unique(z_up) %in% c("HC", "MDD"))) {
    return(factor(z_up, levels = c("HC", "MDD")))
  }
  stop("Unexpected class labels in ", dataset_name, ": ", paste(sort(unique(z)), collapse = ", "))
}
force_numeric_finite <- function(df, genes, dataset_name) {
  for (g in genes) df[[g]] <- suppressWarnings(as.numeric(df[[g]]))
  if (anyNA(df[, genes, drop = FALSE])) {
    stop("Missing or non-numeric gene-expression values detected in ", dataset_name, ".")
  }
  m <- as.matrix(df[, genes, drop = FALSE])
  if (any(!is.finite(m))) stop("Inf/-Inf detected in ", dataset_name, ".")
  df
}
validate_and_align <- function(syn, real, syn_name, expected_syn_rows = 200) {
  if (!"class" %in% names(syn)) stop(syn_name, " must contain a 'class' column.")
  real_genes <- setdiff(names(real), "class")
  syn_genes <- setdiff(names(syn), "class")
  missing_genes <- setdiff(real_genes, syn_genes)
  extra_genes <- setdiff(syn_genes, real_genes)
  if (length(missing_genes) > 0 || length(extra_genes) > 0) {
    stop(
      syn_name, " does not contain exactly the same genes as Real.\n",
      "Missing: ", ifelse(length(missing_genes) == 0, "none", paste(missing_genes, collapse = ", ")),
      "\nExtra: ", ifelse(length(extra_genes) == 0, "none", paste(extra_genes, collapse = ", "))
    )
  }
  if (nrow(syn) != expected_syn_rows) {
    stop(syn_name, " must contain ", expected_syn_rows, " rows; found ", nrow(syn), ".")
  }
  syn <- syn[, c(real_genes, "class"), drop = FALSE]
  syn$class <- fix_class(syn$class, syn_name)
  force_numeric_finite(syn, real_genes, syn_name)
}
rf_oob_accuracy <- function(df, seed_value, trees) {
  m <- ranger(class ~ ., data = df, num.trees = trees, importance = "none", seed = seed_value)
  as.numeric(1 - m$prediction.error)
}
rf_rank <- function(df, seed_value = SEED, trees = RANK_TREES) {
  m <- ranger(
    class ~ ., data = df, num.trees = trees,
    importance = "permutation", scale.permutation.importance = TRUE,
    seed = seed_value
  )
  tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |>
    mutate(rank = row_number())
}
rbo_score <- function(list1, list2, p = 0.95) {
  k <- min(length(list1), length(list2))
  overlap <- numeric(k)
  for (d in seq_len(k)) {
    overlap[d] <- length(intersect(list1[seq_len(d)], list2[seq_len(d)])) / d
  }
  (1 - p) * sum(overlap * p^(0:(k - 1))) + overlap[k] * p^k
}
reorder_within <- function(x, by, within, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by)
}
scale_y_reordered <- function(...) {
  scale_y_discrete(labels = function(x) gsub("___.*$", "", x), ...)
}

## =========================================================
## 4) LOAD + STRICT DATA VALIDATION
## =========================================================
msg("============================================================")
msg("Loading and validating datasets ...")
msg("============================================================")
all_paths <- c(Real = real_path, ARF = arf_path, CTGAN = ctg_path, TVAE = tvae_path)
not_found <- names(all_paths)[!file.exists(all_paths)]
if (length(not_found) > 0) {
  stop("Missing input file(s): ", paste(not_found, collapse = ", "), "\nExpected under: ", base_dir, "\nExpected exact CTGAN filename: synthetic_ctgan.csv")
}
real <- read.csv(real_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"class" %in% names(real)) stop("Real data must contain a 'class' column.")
EXPECTED_REAL_ROWS <- 157
EXPECTED_GENE_FEATURES <- 1783
genes <- setdiff(names(real), "class")
if (nrow(real) != EXPECTED_REAL_ROWS) stop("Real: expected 157 rows, found ", nrow(real), ".")
if (length(genes) != EXPECTED_GENE_FEATURES) stop("Real: expected 1783 genes, found ", length(genes), ".")
real$class <- fix_class(real$class, "Real")
real <- force_numeric_finite(real, genes, "Real")
arf <- validate_and_align(read.csv(arf_path, check.names = FALSE, stringsAsFactors = FALSE), real, "ARF")
ctgan <- validate_and_align(read.csv(ctg_path, check.names = FALSE, stringsAsFactors = FALSE), real, "CTGAN")
tvae <- validate_and_align(read.csv(tvae_path, check.names = FALSE, stringsAsFactors = FALSE), real, "TVAE")
datasets <- list(Real = real, ARF = arf, CTGAN = ctgan, TVAE = tvae)
validation_tbl <- bind_rows(lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]; class_tab <- table(d$class)
  tibble(
    dataset = nm, rows = nrow(d), gene_features = length(genes), columns_total = ncol(d),
    class_HC = unname(class_tab["HC"]), class_MDD = unname(class_tab["MDD"]),
    missing_values = sum(is.na(d)), duplicated_rows = sum(duplicated(d)),
    zero_variance_genes = sum(vapply(d[, genes, drop = FALSE], function(x) sd(x) == 0, logical(1)))
  )
}))
write_csv(validation_tbl, file.path(tables_dir, "00_data_validation.csv"))
print(validation_tbl)
if (any(validation_tbl$zero_variance_genes > 0)) warning("At least one dataset contains zero-variance genes; inspect 00_data_validation.csv.")

## =========================================================
## 5) UTILITY: RF OOB ACCURACY
## =========================================================
msg("\n============================================================")
msg("1/7 Utility: RF OOB accuracy (%d repetitions) ...", UTILITY_REPS)
msg("============================================================")
acc_real <- numeric(UTILITY_REPS)
for (i in seq_len(UTILITY_REPS)) {
  if (i == 1 || i %% 10 == 0 || i == UTILITY_REPS) {
    msg("  Utility Real: repetition %d/%d", i, UTILITY_REPS)
  }
  acc_real[i] <- rf_oob_accuracy(real, SEED + i, UTILITY_TREES)
}
utility_rows <- list(); paired_tests <- list(); gen_names <- c("ARF", "CTGAN", "TVAE")
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; combined <- bind_rows(real, syn)
  msg("Utility: %s ...", g)
  acc_syn <- numeric(UTILITY_REPS)
  acc_comb <- numeric(UTILITY_REPS)
  for (i in seq_len(UTILITY_REPS)) {
    if (i == 1 || i %% 10 == 0 || i == UTILITY_REPS) {
      msg("  %s synthetic-only: repetition %d/%d", g, i, UTILITY_REPS)
    }
    acc_syn[i] <- rf_oob_accuracy(syn, SEED + 10000 * g_idx + i, UTILITY_TREES)
  }
  for (i in seq_len(UTILITY_REPS)) {
    if (i == 1 || i %% 10 == 0 || i == UTILITY_REPS) {
      msg("  Real+%s: repetition %d/%d", g, i, UTILITY_REPS)
    }
    acc_comb[i] <- rf_oob_accuracy(combined, SEED + i, UTILITY_TREES)
  }
  utility_rows[[g]] <- tibble(generator = g, replicate = seq_len(UTILITY_REPS), Real = acc_real, Synthetic = acc_syn, Real_plus_Synthetic = acc_comb)
  tt <- t.test(acc_comb, acc_real, paired = TRUE); diffs <- acc_comb - acc_real
  paired_tests[[g]] <- tibble(
    generator = g, comparison = "Real+Synthetic vs Real", mean_difference = mean(diffs),
    ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], t = unname(tt$statistic), df = unname(tt$parameter),
    p_value = tt$p.value, cohens_dz = unname(tt$statistic) / sqrt(length(diffs))
  )
}
utility_wide <- bind_rows(utility_rows)
utility_long <- utility_wide |>
  pivot_longer(cols = c(Real, Synthetic, Real_plus_Synthetic), names_to = "condition", values_to = "accuracy") |>
  mutate(
    condition = recode(condition, Real = "Real", Synthetic = "Synthetic", Real_plus_Synthetic = "Real+Synthetic"),
    condition = factor(condition, levels = c("Real", "Synthetic", "Real+Synthetic"))
  )
utility_summary <- utility_long |>
  group_by(generator, condition) |>
  summarise(n = n(), mean_accuracy = mean(accuracy), sd_accuracy = sd(accuracy),
            ci_low = mean_accuracy - 1.96 * sd_accuracy / sqrt(n),
            ci_high = mean_accuracy + 1.96 * sd_accuracy / sqrt(n), .groups = "drop")
utility_tests <- bind_rows(paired_tests)
write_csv(utility_wide, file.path(details_dir, "01_utility_replicates.csv"))
write_csv(utility_summary, file.path(tables_dir, "01_utility_summary.csv"))
write_csv(utility_tests, file.path(tables_dir, "01_utility_paired_tests.csv"))
p_utility <- ggplot(utility_long, aes(x = condition, y = accuracy)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
  facet_wrap(~ generator, nrow = 1) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = "RF OOB accuracy", title = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig1_Utility_all_generators.png"), p_utility, width = 11, height = 4.2, dpi = 300)
print(utility_summary)

## =========================================================
## 6) UMAP
## =========================================================
msg("\n============================================================")
msg("2/7 UMAP visualizations ...")
msg("============================================================")
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- UMAP_NEIGHBORS; umap_cfg$min_dist <- UMAP_MIN_DIST
umap_coordinates <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; msg("UMAP: Real vs %s ...", g)
  comb <- bind_rows(mutate(real, source = "Real"), mutate(syn, source = g))
  set.seed(SEED + 5000 + g_idx)
  emb <- as.data.frame(umap(as.matrix(comb[, genes, drop = FALSE]), config = umap_cfg)$layout)
  colnames(emb) <- c("UMAP1", "UMAP2"); emb$source <- comb$source; emb$class <- comb$class; emb$generator <- g
  umap_coordinates[[paste0("Real_vs_", g)]] <- emb
  p_comb <- ggplot(emb, aes(UMAP1, UMAP2, color = source, shape = class)) + geom_point(alpha = 0.75, size = 2) +
    theme_minimal(base_size = 13) + labs(title = paste("Real vs", g), x = "UMAP 1", y = "UMAP 2") + theme(panel.grid.minor = element_blank())
  ggsave(file.path(figures_dir, paste0("Fig2_UMAP_Real_vs_", g, ".png")), p_comb, width = 6, height = 5, dpi = 300)
  set.seed(SEED + 6000 + g_idx)
  emb_syn <- as.data.frame(umap(as.matrix(syn[, genes, drop = FALSE]), config = umap_cfg)$layout)
  colnames(emb_syn) <- c("UMAP1", "UMAP2"); emb_syn$class <- syn$class; emb_syn$source <- g; emb_syn$generator <- g
  umap_coordinates[[paste0(g, "_only")]] <- emb_syn
  p_syn <- ggplot(emb_syn, aes(UMAP1, UMAP2, color = class)) + geom_point(alpha = 0.8, size = 2) +
    theme_minimal(base_size = 13) + labs(title = paste(g, "synthetic"), x = "UMAP 1", y = "UMAP 2") + theme(panel.grid.minor = element_blank())
  ggsave(file.path(figures_dir, paste0("Fig2_UMAP_", g, "_only.png")), p_syn, width = 6, height = 5, dpi = 300)
}
write_csv(bind_rows(umap_coordinates, .id = "panel"), file.path(details_dir, "02_umap_coordinates.csv"))

## =========================================================
## 7) FIDELITY
## =========================================================
msg("\n============================================================")
msg("3/7 Fidelity: KS + Wilcoxon + BH-FDR ...")
msg("============================================================")
fidelity_summaries <- list()
for (g in gen_names) {
  msg("Fidelity: Real vs %s ...", g); syn <- datasets[[g]]
  ks_p <- vapply(genes, function(gene) suppressWarnings(ks.test(real[[gene]], syn[[gene]])$p.value), numeric(1))
  wil_p <- vapply(genes, function(gene) suppressWarnings(wilcox.test(real[[gene]], syn[[gene]])$p.value), numeric(1))
  d <- tibble(generator = g, gene = genes, ks_pvalue = ks_p, wilcoxon_pvalue = wil_p,
              ks_qvalue = p.adjust(ks_p, method = "BH"), wilcoxon_qvalue = p.adjust(wil_p, method = "BH"))
  s <- d |> summarise(generator = first(generator), genes_tested = n(),
                      KS_raw_pct = mean(ks_pvalue < 0.05, na.rm = TRUE) * 100,
                      Wilcoxon_raw_pct = mean(wilcoxon_pvalue < 0.05, na.rm = TRUE) * 100,
                      KS_BH_FDR_pct = mean(ks_qvalue < 0.05, na.rm = TRUE) * 100,
                      Wilcoxon_BH_FDR_pct = mean(wilcoxon_qvalue < 0.05, na.rm = TRUE) * 100)
  fidelity_summaries[[g]] <- s
  write_csv(d, file.path(details_dir, paste0("03_fidelity_gene_level_", g, ".csv")))
}
fidelity_summary <- bind_rows(fidelity_summaries)
write_csv(fidelity_summary, file.path(tables_dir, "03_fidelity_summary.csv"))
fidelity_plot_df <- fidelity_summary |>
  select(generator, KS_raw_pct, Wilcoxon_raw_pct, KS_BH_FDR_pct, Wilcoxon_BH_FDR_pct) |>
  pivot_longer(-generator, names_to = "test", values_to = "percent_significant") |>
  mutate(test = recode(test, KS_raw_pct = "KS raw", Wilcoxon_raw_pct = "Wilcoxon raw", KS_BH_FDR_pct = "KS BH-FDR", Wilcoxon_BH_FDR_pct = "Wilcoxon BH-FDR"))
p_fidelity <- ggplot(fidelity_plot_df, aes(test, percent_significant, fill = generator)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = "Genes significant (%)", fill = "Generator", title = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fidelity_all_generators.png"), p_fidelity, width = 8, height = 5, dpi = 300)
print(fidelity_summary)

## =========================================================
## 8) PRIVACY-RELATED DISTINGUISHABILITY
## =========================================================
msg("\n============================================================")
msg("4/7 Privacy-related distinguishability: %dx%d-fold CV ...", PRIVACY_REPEATS, PRIVACY_FOLDS)
msg("============================================================")
privacy_fold_results <- list(); privacy_summaries <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; msg("Privacy: Real vs %s ...", g)
  dom <- bind_rows(mutate(real, source = "Real"), mutate(syn, source = "Synthetic"))
  dom$source <- factor(dom$source, levels = c("Synthetic", "Real"))
  set.seed(SEED + 7000 + g_idx)
  folds <- vfold_cv(dom, v = PRIVACY_FOLDS, repeats = PRIVACY_REPEATS, strata = source)
  fold_rows <- vector("list", nrow(folds))
  for (i in seq_len(nrow(folds))) {
    train <- analysis(folds$splits[[i]]); test <- assessment(folds$splits[[i]])
    m <- ranger(source ~ . - class, data = train, probability = TRUE, num.trees = PRIVACY_TREES,
                seed = SEED + 8000 + 1000 * g_idx + i)
    pred <- predict(m, data = test)$predictions
    if (!"Real" %in% colnames(pred)) stop("Privacy model for ", g, " did not return probability column 'Real'.")
    prob_real <- pred[, "Real"]
    roc_obj <- pROC::roc(response = test$source, predictor = prob_real,
                         levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)
    auc_i <- as.numeric(pROC::auc(roc_obj))
    fold_rows[[i]] <- tibble(generator = g, fold_id = i, truth = test$source, prob_real = prob_real, fold_auc = auc_i)
  }
  fold_df <- bind_rows(fold_rows)
  pooled_roc <- pROC::roc(response = fold_df$truth, predictor = fold_df$prob_real,
                          levels = c("Synthetic", "Real"), direction = "<", quiet = TRUE)
  pooled_auc <- as.numeric(pROC::auc(pooled_roc)); normalized_auc <- max(pooled_auc, 1 - pooled_auc)
  fold_auc_df <- fold_df |> distinct(generator, fold_id, fold_auc); ci <- ci95(fold_auc_df$fold_auc)
  privacy_fold_results[[g]] <- fold_auc_df
  privacy_summaries[[g]] <- tibble(generator = g, pooled_auc = pooled_auc, normalized_auc = normalized_auc,
                                   mean_fold_auc = mean(fold_auc_df$fold_auc), sd_fold_auc = sd(fold_auc_df$fold_auc),
                                   fold_auc_ci_low = ci[1], fold_auc_ci_high = ci[2], n_folds = nrow(fold_auc_df))
  write_csv(fold_df, file.path(details_dir, paste0("04_privacy_oof_predictions_", g, ".csv")))
}
privacy_folds <- bind_rows(privacy_fold_results); privacy_summary <- bind_rows(privacy_summaries)
write_csv(privacy_folds, file.path(details_dir, "04_privacy_fold_auc.csv"))
write_csv(privacy_summary, file.path(tables_dir, "04_privacy_summary.csv"))
p_privacy <- ggplot(privacy_summary, aes(generator, normalized_auc)) + geom_col(width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + geom_text(aes(label = sprintf("%.3f", normalized_auc)), vjust = -0.5, size = 4) +
  coord_cartesian(ylim = c(0.45, 1.05)) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = "Normalized domain-separator AUC", title = NULL) + theme(panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig3_Privacy_all_generators.png"), p_privacy, width = 6.5, height = 5, dpi = 300)
print(privacy_summary)

## =========================================================
## 9) RF FEATURE-IMPORTANCE ALIGNMENT
## =========================================================
msg("\n============================================================")
msg("5/7 RF feature-importance alignment ...")
msg("============================================================")
rank_real <- rf_rank(real, seed_value = SEED, trees = RANK_TREES)
rf_ranks <- list(Real = rank_real); rf_combined_ranks <- list(); rf_pairwise <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; msg("RF importance: %s ...", g)
  rank_syn <- rf_rank(syn, seed_value = SEED + 9000 + g_idx, trees = RANK_TREES)
  rank_comb <- rf_rank(bind_rows(real, syn), seed_value = SEED, trees = RANK_TREES)
  rf_ranks[[g]] <- rank_syn; rf_combined_ranks[[g]] <- rank_comb
  top_real <- head(rank_real$gene, TOP_K); top_syn <- head(rank_syn$gene, TOP_K); top_comb <- head(rank_comb$gene, TOP_K)
  rf_pairwise[[g]] <- tibble(generator = g,
    comparison = c("Real vs Synthetic", "Real vs Real+Synthetic", "Synthetic vs Real+Synthetic"), TopK = TOP_K,
    n_overlap = c(length(intersect(top_real, top_syn)), length(intersect(top_real, top_comb)), length(intersect(top_syn, top_comb))),
    jaccard = c(jaccard(top_real, top_syn), jaccard(top_real, top_comb), jaccard(top_syn, top_comb)))
  write_csv(rank_syn, file.path(details_dir, paste0("05_rf_importance_", g, ".csv")))
  write_csv(rank_comb, file.path(details_dir, paste0("05_rf_importance_RealPlus", g, ".csv")))
}
write_csv(rank_real, file.path(details_dir, "05_rf_importance_Real.csv"))
rf_alignment_summary <- bind_rows(rf_pairwise)
write_csv(rf_alignment_summary, file.path(tables_dir, "05_rf_feature_alignment_summary.csv"))
rf_top20_plot_df <- bind_rows(lapply(names(rf_ranks), function(nm) rf_ranks[[nm]] |> slice_head(n = TOP_BAR) |> mutate(dataset = nm))) |>
  mutate(gene_facet = reorder_within(gene, importance, dataset))
p_rf <- ggplot(rf_top20_plot_df, aes(importance, gene_facet)) + geom_col() + facet_wrap(~ dataset, scales = "free_y", ncol = 2) +
  scale_y_reordered() + theme_minimal(base_size = 11) + labs(x = "Permutation importance", y = NULL, title = NULL) + theme(panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig4_RF_importance_top20_all_generators.png"), p_rf, width = 11, height = 10, dpi = 300)
print(rf_alignment_summary)

## =========================================================
## 10) NPDR FEATURE ALIGNMENT
## =========================================================
msg("\n============================================================")
msg("6/7 NPDR feature alignment ...")
msg("============================================================")
run_npdr_rank <- function(df, dataset_label) {
  msg("NPDR: %s ...", dataset_label)
  X <- as.data.frame(df[, genes, drop = FALSE]); y <- droplevels(df$class)
  res <- npdr::npdr(outcome = y, dataset = X, regression.type = "binomial", attr.diff.type = "numeric-abs",
                    nbd.method = NPDR_METHOD, nbd.metric = NPDR_METRIC, knn = NPDR_KNN,
                    padj.method = NPDR_PADJ, verbose = FALSE)
  required_cols <- c("att", "pval.adj", "pval.att", "beta.raw.att")
  missing_cols <- setdiff(required_cols, names(res))
  if (length(missing_cols) > 0) {
    stop("Unexpected NPDR output for ", dataset_label, ". Missing: ", paste(missing_cols, collapse = ", "),
         "\nObserved columns: ", paste(names(res), collapse = ", "))
  }
  ranked <- as_tibble(res) |>
    transmute(gene = as.character(att), pval_adj = as.numeric(pval.adj), pval_raw = as.numeric(pval.att),
              beta_raw = as.numeric(beta.raw.att), beta_z = if ("beta.Z.att" %in% names(res)) as.numeric(beta.Z.att) else NA_real_) |>
    arrange(pval_adj, pval_raw, desc(beta_raw), gene) |>
    mutate(rank = row_number())
  if (!setequal(ranked$gene, genes)) stop("NPDR gene names do not match canonical genes for ", dataset_label, ".")
  ranked
}
npdr_real <- run_npdr_rank(real, "Real")
write_csv(npdr_real, file.path(details_dir, "06_npdr_importance_Real.csv"))
npdr_pairwise <- list()
for (g in gen_names) {
  syn <- datasets[[g]]; comb <- bind_rows(real, syn)
  npdr_syn <- run_npdr_rank(syn, g); npdr_comb <- run_npdr_rank(comb, paste0("Real+", g))
  write_csv(npdr_syn, file.path(details_dir, paste0("06_npdr_importance_", g, ".csv")))
  write_csv(npdr_comb, file.path(details_dir, paste0("06_npdr_importance_RealPlus", g, ".csv")))
  top_real <- head(npdr_real$gene, TOP_K); top_syn <- head(npdr_syn$gene, TOP_K); top_comb <- head(npdr_comb$gene, TOP_K)
  npdr_pairwise[[g]] <- tibble(generator = g,
    comparison = c("Real vs Synthetic", "Real vs Real+Synthetic", "Synthetic vs Real+Synthetic"), TopK = TOP_K,
    n_overlap = c(length(intersect(top_real, top_syn)), length(intersect(top_real, top_comb)), length(intersect(top_syn, top_comb))),
    jaccard = c(jaccard(top_real, top_syn), jaccard(top_real, top_comb), jaccard(top_syn, top_comb)))
}
npdr_alignment_summary <- bind_rows(npdr_pairwise)
write_csv(npdr_alignment_summary, file.path(tables_dir, "06_npdr_feature_alignment_summary.csv"))
p_npdr <- ggplot(npdr_alignment_summary, aes(comparison, jaccard, fill = generator)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = paste0("NPDR Top-", TOP_K, " Jaccard"), fill = "Generator", title = NULL) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig5_NPDR_alignment_all_generators.png"), p_npdr, width = 9, height = 5, dpi = 300)
print(npdr_alignment_summary)

## =========================================================
## 11) RANK AGREEMENT + DILUTION
## =========================================================
msg("\n============================================================")
msg("Ranking agreement + dilution ...")
msg("============================================================")
calculate_rank_metrics <- function(real_rank, augmented_rank, comparison_name, K = TOP_K, rbo_p = 0.95) {
  merged <- real_rank |> select(gene, rank_real = rank) |>
    inner_join(augmented_rank |> select(gene, rank_aug = rank), by = "gene")
  if (nrow(merged) != length(genes)) stop("Ranking merge failed for ", comparison_name, ".")
  spearman_rho <- cor(merged$rank_real, merged$rank_aug, method = "spearman")
  kendall_tau <- cor(merged$rank_real, merged$rank_aug, method = "kendall")
  real_list <- real_rank$gene; aug_list <- augmented_rank$gene
  rbo <- rbo_score(real_list, aug_list, p = rbo_p)
  top_real <- head(real_list, K); top_aug <- head(aug_list, K); jac <- jaccard(top_real, top_aug)
  detail <- merged |> filter(gene %in% top_real) |>
    mutate(signed_rank_shift = rank_aug - rank_real, absolute_rank_shift = abs(signed_rank_shift), dropped = rank_aug > K) |>
    arrange(rank_real)
  summary <- tibble(Comparison = comparison_name, Spearman_rho = spearman_rho, Kendall_tau = kendall_tau,
                    RBO_p_0.95 = rbo, Jaccard_K50 = jac,
                    Mean_absolute_rank_shift = mean(detail$absolute_rank_shift),
                    Median_absolute_rank_shift = median(detail$absolute_rank_shift),
                    Top50_dropout_pct = mean(detail$dropped) * 100,
                    Mean_signed_rank_shift = mean(detail$signed_rank_shift),
                    Median_signed_rank_shift = median(detail$signed_rank_shift))
  list(summary = summary, detail = detail)
}
rank_metric_results <- list()
for (g in gen_names) {
  res <- calculate_rank_metrics(rank_real, rf_combined_ranks[[g]], paste0("Real vs Real+", g), K = TOP_K, rbo_p = 0.95)
  rank_metric_results[[g]] <- res
  write_csv(res$detail, file.path(details_dir, paste0("07_dilution_detail_", g, "_K", TOP_K, ".csv")))
}
ranking_dilution_summary <- bind_rows(lapply(rank_metric_results, `[[`, "summary"))
write_csv(ranking_dilution_summary, file.path(tables_dir, "07_ranking_dilution_summary_K50.csv"))
dilution_plot_df <- ranking_dilution_summary |>
  transmute(generator = gsub("^Real vs Real\\+", "", Comparison),
            `Mean absolute rank shift` = Mean_absolute_rank_shift,
            `Median absolute rank shift` = Median_absolute_rank_shift,
            `Top-50 dropout (%)` = Top50_dropout_pct) |>
  pivot_longer(-generator, names_to = "metric", values_to = "value")
p_dilution <- ggplot(dilution_plot_df, aes(generator, value)) + geom_col(width = 0.65) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = NULL, title = NULL) + theme(panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig6_Dilution_all_generators.png"), p_dilution, width = 11, height = 4.3, dpi = 300)
print(ranking_dilution_summary)

## =========================================================
## 12) BOOTSTRAP FEATURE-SELECTION STABILITY
## =========================================================
msg("\n============================================================")
msg("7/7 Bootstrap RF Top-%d stability (%d replicates) ...", TOP_K, STABILITY_B)
msg("This is usually the slowest section.")
msg("============================================================")
get_topk_rf <- function(df, seed_value, k = TOP_K) {
  m <- ranger(class ~ ., data = df, num.trees = STABILITY_TREES,
              importance = "permutation", scale.permutation.importance = TRUE, seed = seed_value)
  tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |> slice_head(n = k) |> pull(gene)
}
bootstrap_self_stability <- function(df, label, B = STABILITY_B, seed_offset = 0) {
  n <- nrow(df); out <- numeric(B)
  for (b in seq_len(B)) {
    if (b %% 10 == 0 || b == 1 || b == B) msg("  %s self stability: %d/%d", label, b, B)
    set.seed(SEED + seed_offset + 100 * b + 1); idx1 <- sample.int(n, size = n, replace = TRUE)
    set.seed(SEED + seed_offset + 100 * b + 2); idx2 <- sample.int(n, size = n, replace = TRUE)
    top1 <- get_topk_rf(df[idx1, , drop = FALSE], SEED + seed_offset + 100 * b + 3)
    top2 <- get_topk_rf(df[idx2, , drop = FALSE], SEED + seed_offset + 100 * b + 4)
    out[b] <- jaccard(top1, top2)
  }
  tibble(comparison = label, replicate = seq_len(B), jaccard = out)
}
bootstrap_cross_stability <- function(real_df, syn_df, label, B = STABILITY_B, seed_offset = 0) {
  nr <- nrow(real_df); ns <- nrow(syn_df); out <- numeric(B)
  for (b in seq_len(B)) {
    if (b %% 10 == 0 || b == 1 || b == B) msg("  %s: %d/%d", label, b, B)
    set.seed(SEED + seed_offset + 100 * b + 1); idx_r <- sample.int(nr, size = nr, replace = TRUE)
    set.seed(SEED + seed_offset + 100 * b + 2); idx_s <- sample.int(ns, size = ns, replace = TRUE)
    top_r <- get_topk_rf(real_df[idx_r, , drop = FALSE], SEED + seed_offset + 100 * b + 3)
    top_s <- get_topk_rf(syn_df[idx_s, , drop = FALSE], SEED + seed_offset + 100 * b + 4)
    out[b] <- jaccard(top_r, top_s)
  }
  tibble(comparison = label, replicate = seq_len(B), jaccard = out)
}
stability_parts <- list()
stability_parts[["Real"]] <- bootstrap_self_stability(real, "Real", seed_offset = 100000)
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]
  stability_parts[[g]] <- bootstrap_self_stability(datasets[[g]], g, seed_offset = 200000 + 10000 * g_idx)
}
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]
  stability_parts[[paste0("Real_vs_", g)]] <- bootstrap_cross_stability(real, datasets[[g]], paste0("Real vs ", g), seed_offset = 300000 + 10000 * g_idx)
}
stability_replicates <- bind_rows(stability_parts)
stability_summary <- stability_replicates |>
  group_by(comparison) |>
  summarise(B = n(), mean_jaccard = mean(jaccard), sd_jaccard = sd(jaccard),
            ci_low = mean_jaccard - 1.96 * sd_jaccard / sqrt(B),
            ci_high = mean_jaccard + 1.96 * sd_jaccard / sqrt(B), .groups = "drop")
write_csv(stability_replicates, file.path(details_dir, "08_stability_replicates.csv"))
write_csv(stability_summary, file.path(tables_dir, "08_stability_summary.csv"))
p_stability <- ggplot(stability_summary, aes(comparison, mean_jaccard)) + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18) + theme_minimal(base_size = 13) +
  labs(x = NULL, y = paste0("Bootstrap Top-", TOP_K, " Jaccard"), title = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "Fig7_Stability_all_generators.png"), p_stability, width = 9, height = 5.5, dpi = 300)
print(stability_summary)

## =========================================================
## 13) MASTER SUMMARY FOR REVIEWER COMMENT 2
## =========================================================
msg("\n============================================================")
msg("Creating compact master summary for Reviewer 1, Comment 2 ...")
msg("============================================================")
master_summary <- tibble(generator = gen_names) |>
  left_join(utility_summary |> filter(condition == "Synthetic") |>
              transmute(generator, utility_synthetic_mean = mean_accuracy, utility_synthetic_sd = sd_accuracy), by = "generator") |>
  left_join(utility_summary |> filter(condition == "Real+Synthetic") |>
              transmute(generator, utility_realplus_mean = mean_accuracy, utility_realplus_sd = sd_accuracy), by = "generator") |>
  left_join(fidelity_summary |> select(generator, KS_BH_FDR_pct, Wilcoxon_BH_FDR_pct), by = "generator") |>
  left_join(privacy_summary |> select(generator, privacy_pooled_auc = pooled_auc, privacy_normalized_auc = normalized_auc), by = "generator") |>
  left_join(rf_alignment_summary |> filter(comparison == "Real vs Synthetic") |>
              transmute(generator, RF_Top50_overlap = n_overlap, RF_Top50_Jaccard = jaccard), by = "generator") |>
  left_join(npdr_alignment_summary |> filter(comparison == "Real vs Synthetic") |>
              transmute(generator, NPDR_Top50_overlap = n_overlap, NPDR_Top50_Jaccard = jaccard), by = "generator") |>
  left_join(ranking_dilution_summary |>
              transmute(generator = gsub("^Real vs Real\\+", "", Comparison), Spearman_rho, Kendall_tau, RBO_p_0.95,
                        Jaccard_K50, Mean_absolute_rank_shift, Median_absolute_rank_shift, Top50_dropout_pct), by = "generator") |>
  left_join(stability_summary |> filter(comparison %in% gen_names) |>
              transmute(generator = comparison, self_stability_mean_jaccard = mean_jaccard), by = "generator") |>
  left_join(stability_summary |> filter(grepl("^Real vs ", comparison)) |>
              transmute(generator = sub("^Real vs ", "", comparison),
                        real_vs_synthetic_stability_mean_jaccard = mean_jaccard), by = "generator")
write_csv(master_summary, file.path(tables_dir, "09_MASTER_SUMMARY_Reviewer1_Comment2.csv"))
print(master_summary)

## =========================================================
## 14) REPRODUCIBILITY LOGS
## =========================================================
settings_tbl <- tibble(
  setting = c("QUICK_TEST", "SEED", "TOP_K", "UTILITY_REPS", "UTILITY_TREES", "RANK_TREES",
              "PRIVACY_FOLDS", "PRIVACY_REPEATS", "PRIVACY_TREES", "STABILITY_B", "STABILITY_TREES",
              "NPDR_METHOD", "NPDR_METRIC", "NPDR_KNN", "NPDR_PADJ", "UMAP_NEIGHBORS", "UMAP_MIN_DIST"),
  value = as.character(c(QUICK_TEST, SEED, TOP_K, UTILITY_REPS, UTILITY_TREES, RANK_TREES,
                         PRIVACY_FOLDS, PRIVACY_REPEATS, PRIVACY_TREES, STABILITY_B, STABILITY_TREES,
                         NPDR_METHOD, NPDR_METRIC, NPDR_KNN, NPDR_PADJ, UMAP_NEIGHBORS, UMAP_MIN_DIST))
)
write_csv(settings_tbl, file.path(logs_dir, "run_settings.csv"))
capture.output(sessionInfo(), file = file.path(logs_dir, "sessionInfo.txt"))
writeLines(c(
  "Unified ARF / CTGAN / TVAE evaluation", "",
  "Key manuscript-ready tables:",
  "  tables/01_utility_summary.csv",
  "  tables/03_fidelity_summary.csv",
  "  tables/04_privacy_summary.csv",
  "  tables/05_rf_feature_alignment_summary.csv",
  "  tables/06_npdr_feature_alignment_summary.csv",
  "  tables/07_ranking_dilution_summary_K50.csv",
  "  tables/08_stability_summary.csv", "",
  "Most useful compact file:",
  "  tables/09_MASTER_SUMMARY_Reviewer1_Comment2.csv", "",
  paste0("QUICK_TEST = ", QUICK_TEST),
  paste0("Run completed: ", Sys.time())
), file.path(out_dir, "README_RESULTS.txt"))

## =========================================================
## 15) DONE
## =========================================================
msg("\n============================================================")
msg("ALL ANALYSES COMPLETED SUCCESSFULLY.")
msg("============================================================")
msg("Results folder:")
msg("%s", out_dir)
msg("\nFIRST FILE TO SEND BACK FOR INTERPRETATION:")
msg("%s", file.path(tables_dir, "09_MASTER_SUMMARY_Reviewer1_Comment2.csv"))
msg("============================================================")
