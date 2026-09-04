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

## Shared color palette for ARF / CTGAN / TVAE, used consistently in every
## figure in the manuscript. ARF/CTGAN match the colors already used in the
## published figures; TVAE uses a bold, colorblind-safe amber/gold so it
## reads clearly next to red and teal in print and in grayscale.
GENERATOR_COLORS <- c(
  ARF   = "#F8766D",  # coral red  (same as current manuscript figures)
  CTGAN = "#00BFC4",  # teal       (same as current manuscript figures)
  TVAE  = "#E5A500"   # bold amber/gold
)

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
required_cran <- c("tidyverse", "ranger", "umap", "pROC", "rsample", "ragg", "patchwork", "cowplot")
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
  library(ragg)
})

## Shared theme: Times New Roman, horizontal axis labels, used in every figure.
FONT_FAMILY <- "Times New Roman"
THEME_BASE <- theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    strip.text = element_text(family = FONT_FAMILY, face = "bold", size = 12),
    legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12),
    legend.title = element_text(family = FONT_FAMILY, face = "bold", size = 12)
  )
save_fig <- function(filename, plot, width, height) {
  ggsave(filename, plot, width = width, height = height, dpi = 600, device = ragg::agg_png, bg = "white")
}
set.seed(SEED)

## =========================================================
## 3) HELPERS
## =========================================================
msg <- function(...) cat(sprintf(...), "\n")
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
  stop("Missing input file(s): ", paste(not_found, collapse = ", "), "\nExpected under: ", base_dir)
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
acc_real <- vapply(seq_len(UTILITY_REPS), function(i) rf_oob_accuracy(real, SEED + i, UTILITY_TREES), numeric(1))
utility_rows <- list(); paired_tests <- list(); gen_names <- c("ARF", "CTGAN", "TVAE")
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; combined <- bind_rows(real, syn)
  msg("Utility: %s ...", g)
  acc_syn <- vapply(seq_len(UTILITY_REPS), function(i) rf_oob_accuracy(syn, SEED + 10000 * g_idx + i, UTILITY_TREES), numeric(1))
  acc_comb <- vapply(seq_len(UTILITY_REPS), function(i) rf_oob_accuracy(combined, SEED + i, UTILITY_TREES), numeric(1))
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
p_utility <- ggplot(utility_long, aes(x = condition, y = accuracy, fill = generator)) +
  geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 2.5, stroke = 1, color = "black") +
  facet_wrap(~ generator, nrow = 1) +
  scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
  THEME_BASE +
  labs(x = NULL, y = "RF OOB accuracy", title = NULL) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    legend.position = "top",
    legend.justification = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.margin = margin(3, 8, 3, 8)
  )
save_fig(file.path(figures_dir, "Fig1_Utility_all_generators.png"), p_utility, width = 13, height = 6.2)
print(utility_summary)

## =========================================================
## 6) UMAP
## =========================================================
msg("\n============================================================")
msg("2/7 UMAP visualizations ...")
msg("============================================================")
library(patchwork)
umap_cfg <- umap.defaults; umap_cfg$n_neighbors <- UMAP_NEIGHBORS; umap_cfg$min_dist <- UMAP_MIN_DIST
umap_coordinates <- list()
for (g_idx in seq_along(gen_names)) {
  g <- gen_names[g_idx]; syn <- datasets[[g]]; msg("UMAP: Real vs %s ...", g)
  comb <- bind_rows(mutate(real, source = "Real"), mutate(syn, source = g))
  set.seed(SEED + 5000 + g_idx)
  emb <- as.data.frame(umap(as.matrix(comb[, genes, drop = FALSE]), config = umap_cfg)$layout)
  colnames(emb) <- c("UMAP1", "UMAP2"); emb$source <- comb$source; emb$class <- comb$class; emb$generator <- g
  umap_coordinates[[paste0("Real_vs_", g)]] <- emb
  set.seed(SEED + 6000 + g_idx)
  emb_syn <- as.data.frame(umap(as.matrix(syn[, genes, drop = FALSE]), config = umap_cfg)$layout)
  colnames(emb_syn) <- c("UMAP1", "UMAP2"); emb_syn$class <- syn$class; emb_syn$source <- g; emb_syn$generator <- g
  umap_coordinates[[paste0(g, "_only")]] <- emb_syn
}
write_csv(bind_rows(umap_coordinates, .id = "panel"), file.path(details_dir, "02_umap_coordinates.csv"))

## Single combined 2x3 figure (a-f), one shared top-right legend.
## NOTE: patchwork's guides="collect" can fail to merge identical per-panel
## legends into one (it sometimes duplicates them instead), so the legend is
## built once from a single reference panel and manually placed with cowplot.
library(cowplot)
ALL_UMAP_COLORS <- c(GENERATOR_COLORS, Real = "grey40")
make_umap_panel <- function(df, show_legend = FALSE, extra_top = 0, extra_bottom = 0) {
  df$source <- factor(df$source, levels = c("Real", "ARF", "CTGAN", "TVAE"))
  df$class <- factor(df$class, levels = c("HC", "MDD"))
  ggplot(df, aes(UMAP1, UMAP2, color = source, shape = class)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(values = ALL_UMAP_COLORS, name = NULL, drop = FALSE, breaks = c("ARF", "CTGAN", "TVAE")) +
    scale_shape_manual(values = c(HC = 17, MDD = 15), name = NULL, drop = FALSE) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 15, size = 4)),
      shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3))
    ) +
    THEME_BASE + labs(x = "V1", y = "V2") +
    theme(panel.grid.minor = element_blank(),
          legend.box = "vertical",
          legend.position = if (show_legend) "top" else "none",
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
          legend.margin = margin(8, 14, 8, 14),
          legend.key.size = unit(0.65, "cm"),
          legend.spacing.x = unit(0.3, "cm"),
          legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12, margin = margin(r = 10)),
          plot.margin = margin(t = 5 + extra_top, b = 5 + extra_bottom, l = 5, r = 5))
}
top_row_data <- lapply(gen_names, function(g) umap_coordinates[[paste0("Real_vs_", g)]])
bottom_row_data <- lapply(gen_names, function(g) umap_coordinates[[paste0(g, "_only")]] |> mutate(source = g))

legend_dummy <- expand.grid(
  UMAP1 = 0, UMAP2 = 0,
  source = c("ARF", "CTGAN", "TVAE"),
  class = c("HC", "MDD")
)
shared_legend <- cowplot::get_legend(make_umap_panel(legend_dummy, show_legend = TRUE))

## Arrange as 3 rows x 2 cols, grouped by generator, matching the original
## Fig.2 layout: (a) Real+ARF, (b) ARF only, (c) Real+CTGAN, (d) CTGAN only,
## (e) Real+TVAE, (f) TVAE only.
panel_rows <- lapply(seq_along(gen_names), function(i) {
  left  <- make_umap_panel(top_row_data[[i]],    extra_bottom = 8)
  right <- make_umap_panel(bottom_row_data[[i]], extra_bottom = 8)
  left | right
})

panels_grid <- (panel_rows[[1]] / panel_rows[[2]] / panel_rows[[3]]) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = 12),
        plot.tag.position = "top")

legend_row <- cowplot::plot_grid(NULL, shared_legend, ncol = 2, rel_widths = c(0.55, 0.45))
p_umap_combined <- cowplot::plot_grid(legend_row, panels_grid, ncol = 1, rel_heights = c(0.11, 1))
save_fig(file.path(figures_dir, "Fig2_UMAP_all_generators.png"), p_umap_combined, width = 12, height = 13.5)
## (a) Real+ARF (b) ARF only (c) Real+CTGAN (d) CTGAN only (e) Real+TVAE (f) TVAE only


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
privacy_plot_df <- privacy_summary |>
  mutate(generator = factor(generator, levels = c("ARF", "CTGAN", "TVAE")))
p_privacy <- ggplot(privacy_plot_df, aes(generator, normalized_auc, fill = generator)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.3f", normalized_auc)),
            vjust = -0.5, size = 4, family = FONT_FAMILY, fontface = "bold") +
  scale_fill_manual(values = GENERATOR_COLORS, guide = "none") +
  coord_cartesian(ylim = c(0, 1.08)) + THEME_BASE +
  labs(x = NULL, y = "Normalized domain-separator AUC", title = NULL) +
  theme(panel.grid.minor = element_blank())
save_fig(file.path(figures_dir, "Fig3_Privacy_all_generators.png"), p_privacy, width = 6.5, height = 5.5)
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
rf_datasets <- c("Real", "ARF", "CTGAN", "TVAE")
RF_PLOT_COLORS <- c(GENERATOR_COLORS, Real = "grey40")
rf_data <- lapply(rf_datasets, function(nm) rf_ranks[[nm]] |> slice_head(n = TOP_BAR) |> mutate(dataset = nm))
names(rf_data) <- rf_datasets

make_rf_panel <- function(df, show_legend = FALSE) {
  df <- df |> mutate(gene_facet = reorder(gene, importance),
                      dataset = factor(dataset, levels = rf_datasets))
  ggplot(df, aes(importance, gene_facet, fill = dataset)) +
    geom_col() +
    scale_fill_manual(values = RF_PLOT_COLORS, name = NULL, drop = FALSE) +
    THEME_BASE + labs(x = "Permutation importance", y = "Gene") +
    theme(panel.grid.minor = element_blank(),
          legend.position = if (show_legend) "top" else "none",
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
          legend.margin = margin(8, 18, 8, 18),
          legend.key.size = unit(0.65, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12, margin = margin(r = 12)))
}

## Shared legend built from a dummy dataset containing all 4 categories,
## so every color swatch renders correctly (see the Fig.2 legend fix).
legend_dummy <- tibble(gene = paste0("g", seq_along(rf_datasets)), importance = 1, dataset = rf_datasets)
shared_legend <- cowplot::get_legend(make_rf_panel(legend_dummy, show_legend = TRUE))

rf_panels <- lapply(rf_datasets, function(nm) make_rf_panel(rf_data[[nm]]))
rf_grid <- (rf_panels[[1]] | rf_panels[[2]]) / (rf_panels[[3]] | rf_panels[[4]]) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = 12),
        plot.tag.position = "top")
legend_row <- cowplot::plot_grid(NULL, shared_legend, ncol = 2, rel_widths = c(0.35, 0.65))
p_rf <- cowplot::plot_grid(legend_row, rf_grid, ncol = 1, rel_heights = c(0.08, 1))
save_fig(file.path(figures_dir, "Fig4_RF_importance_top20_all_generators.png"), p_rf, width = 13, height = 12)
## (a) Real (b) ARF (c) CTGAN (d) TVAE
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
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + THEME_BASE +
  scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
  labs(x = NULL, y = paste0("NPDR Top-", TOP_K, " Jaccard"), fill = "Generator", title = NULL) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "right")
save_fig(file.path(figures_dir, "Fig5_NPDR_alignment_all_generators.png"), p_npdr, width = 9, height = 5)
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
  pivot_longer(-generator, names_to = "metric", values_to = "value") |>
  mutate(generator = factor(generator, levels = c("ARF", "CTGAN", "TVAE")),
         metric = factor(metric, levels = c("Mean absolute rank shift", "Median absolute rank shift", "Top-50 dropout (%)")))
p_dilution <- ggplot(dilution_plot_df, aes(generator, value, fill = generator)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.5, size = 4, family = FONT_FAMILY, fontface = "bold") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
  THEME_BASE +
  labs(x = NULL, y = "Value", title = NULL) +
  theme(panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top", legend.justification = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.margin = margin(3, 8, 3, 8))
save_fig(file.path(figures_dir, "Fig6_Dilution_all_generators.png"), p_dilution, width = 12, height = 5)
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
STABILITY_COLORS <- c(GENERATOR_COLORS, Real = "grey50")
STABILITY_ORDER <- c("Real", "ARF", "Real vs ARF", "CTGAN", "Real vs CTGAN", "TVAE", "Real vs TVAE")
stability_plot_df <- stability_replicates |>
  mutate(fill_group = case_when(
    comparison == "Real" ~ "Real",
    comparison %in% gen_names ~ comparison,
    TRUE ~ sub("^Real vs ", "", comparison)
  )) |>
  mutate(comparison = factor(comparison, levels = STABILITY_ORDER),
         fill_group = factor(fill_group, levels = c("Real", gen_names)))
p_stability <- ggplot(stability_plot_df, aes(comparison, jaccard, fill = fill_group)) +
  geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = STABILITY_COLORS, name = NULL, breaks = gen_names) +
  THEME_BASE +
  labs(x = NULL, y = paste0("Bootstrap Top-", TOP_K, " Jaccard"), title = NULL) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        legend.position = "top", legend.justification = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.margin = margin(3, 8, 3, 8))
save_fig(file.path(figures_dir, "Fig7_Stability_all_generators.png"), p_stability, width = 12, height = 6)
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
