############################################################
# 06c_sim_evaluate.R   (known-signal benchmark, step 3)
#
# Question. In simulated data with KNOWN informative genes, does the dilution of
# the reference Top-50 (Table 3 summaries) track the loss of true signal after
# augmentation, and does it add information beyond rank-agreement / overlap
# measures (Spearman, RBO, Top-K Jaccard)?
#
# For each simulated data set ("real" = sim_XX.csv, 157 subjects x 1,800 genes):
#   Reference   random forest (permutation importance, 3000 trees) on the simulated data
#   Observed    same random forest on simulated data + 200 ARF / CTGAN / TVAE subjects
#   Control A   simulated data + 200 rows resampled from the simulated data (class counts of
#               the matching synthetic set): a harmless perturbation
#   Control B   simulated data re-ranked with other RF seeds (pure ranking noise)
#   Control C   simulated data + 200 resampled rows with PERMUTED class labels: a
#               label-permuted perturbation (positive control; the original rows are unchanged)
# Every ranking is scored (i) against the reference ranking with the manuscript metrics
# (Spearman, Kendall, RBO, Top-K Jaccard, mean/median absolute rank shift, dropout of the
# reference Top-50) and (ii) against the TRUTH: number of true main-effect genes in the
# Top-50, loss relative to the reference, true-signal dropout, noise genes entering the
# Top-50, and the AUROC of the importance scores for true main-effect vs noise genes.
# Interaction-only signals are reported separately (random-forest permutation importance
# rarely recovers genes that act only through interactions).
#
# Usage:  Rscript 06c_sim_evaluate.R        (KS_QUICK=1 for a smoke test)
# Output: results/ (per_ranking_metrics.csv, table_condition_summary.csv,
#         table_paired_vs_controls.csv, table_metric_vs_truth_loss.csv,
#         table_incremental_R2.csv, table_detection_auc.csv, figures, run_info.txt)
############################################################
ks_dir <- Sys.getenv("KS_DIR", unset = ".")
QUICK <- identical(Sys.getenv("KS_QUICK", "0"), "1")
KS_SIMS <- Sys.getenv("KS_SIMS", "")          # e.g. KS_SIMS=1,2 evaluates only these simulations (results_PARTIAL)
SFX <- if (QUICK) "_QUICK" else if (nzchar(KS_SIMS)) "_PARTIAL" else ""
draws_dir <- file.path(ks_dir, paste0("draws", if (QUICK) "_QUICK" else "")); res_dir <- file.path(ks_dir, paste0("results", SFX))
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
GENS <- c("ARF", "CTGAN", "TVAE"); TOP_K <- 50L; N_BOOT <- 2000L
RANK_TREES <- if (QUICK) 300L else 3000L
suppressPackageStartupMessages({ library(ranger); library(dplyr); library(tidyr) })
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
seeds <- read.csv(file.path(ks_dir, "seeds.csv"), stringsAsFactors = FALSE)
truth <- read.csv(file.path(ks_dir, "truth.csv"), stringsAsFactors = FALSE)
get_seed <- function(i, role) { s <- seeds$seed[seeds$sim == i & seeds$role == role]; stopifnot(length(s) == 1); s }
N_CTRL_B <- sum(grepl("^ctrlB_rank_", unique(seeds$role))); N_CTRL_C <- sum(grepl("^ctrlC_rank_", unique(seeds$role)))
if (QUICK) N_CTRL_B <- 3L
sims <- sort(unique(seeds$sim))
sims <- sims[file.exists(file.path(ks_dir, "sims", sprintf("sim_%02d.csv", sims)))]
if (QUICK) sims <- sims[1:min(2, length(sims))]
if (nzchar(KS_SIMS)) sims <- intersect(sims, as.integer(strsplit(KS_SIMS, ",")[[1]]))

rf_rank <- function(df, seed) {
  m <- ranger(class ~ ., data = df, num.trees = RANK_TREES, importance = "permutation",
              scale.permutation.importance = TRUE, seed = seed)
  tibble::tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |> mutate(rank = row_number())
}
# --- metrics identical to 05c_evaluate_repeated_draws.R ---
rbo_fast <- function(ra, rb, p) {
  k <- length(ra); X <- cumsum(tabulate(pmax(ra, rb), nbins = k)); d <- seq_len(k); agree <- X / d
  (1 - p) * sum(agree * p^(d - 1)) + agree[k] * p^k
}
rank_metrics <- function(ref, aug, K = TOP_K) {
  ra <- ref$rank; rb <- aug$rank[match(ref$gene, aug$gene)]; stopifnot(!anyNA(rb))
  jac <- function(k) { a <- ref$gene[ra <= k]; b <- ref$gene[rb <= k]; length(intersect(a, b)) / length(union(a, b)) }
  top <- ra <= K; sh <- rb[top] - ra[top]
  tibble::tibble(spearman = cor(ra, rb, method = "spearman"), kendall = cor(ra, rb, method = "kendall"),
    `rbo_p0.95` = rbo_fast(ra, rb, 0.95), `rbo_p0.98` = rbo_fast(ra, rb, 0.98),
    jaccard_top25 = jac(25), jaccard_top50 = jac(50), jaccard_top100 = jac(100),
    mean_abs_shift = mean(abs(sh)), median_abs_shift = median(abs(sh)), dropout_pct = 100 * mean(rb[top] > K))
}
# --- truth-based metrics ---
auc_pair <- function(x1, x0) mean(outer(x1, x0, ">")) + 0.5 * mean(outer(x1, x0, "=="))
truth_metrics <- function(rk, ref, g_main, g_int, g_noise) {
  top <- rk$gene[rk$rank <= TOP_K]; topref <- ref$gene[ref$rank <= TOP_K]
  ref_true_main <- intersect(topref, g_main)
  tibble::tibble(
    tp50_main = sum(top %in% g_main), tp50_int = sum(top %in% g_int), noise50 = sum(top %in% g_noise),
    true_main_dropout = if (length(ref_true_main)) mean(!(ref_true_main %in% top)) else NA_real_,
    noise_entries = sum(top %in% g_noise & !(top %in% topref)),
    auroc_main = auc_pair(rk$importance[rk$gene %in% g_main], rk$importance[rk$gene %in% g_noise]))
}
read_syn <- function(g, i) {
  f <- file.path(draws_dir, g, sprintf("sim_%02d.csv", i)); if (!file.exists(f)) return(NULL)
  s <- read.csv(f, check.names = FALSE); s$class <- factor(as.character(s$class), levels = c("0", "1")); s
}

rows <- list(); k <- 0L
for (i in sims) {
  real <- read.csv(file.path(ks_dir, "sims", sprintf("sim_%02d.csv", i)), check.names = FALSE)
  genes <- setdiff(names(real), "class"); stopifnot(nrow(real) == 157, setequal(unique(real$class), c(0, 1)))
  real$class <- factor(as.character(real$class), levels = c("0", "1"))
  tr <- truth[truth$sim == i, ]; g_main <- tr$gene[tr$type == "main"]; g_int <- tr$gene[tr$type == "interaction"]
  g_noise <- setdiff(genes, c(g_main, g_int))
  stopifnot(all(c(g_main, g_int) %in% genes), length(g_main) == 45, length(g_int) == 45)
  ref <- rf_rank(real, get_seed(i, "rank_ref"))
  tm_ref <- truth_metrics(ref, ref, g_main, g_int, g_noise)
  add <- function(cond, gen, rep, rk, nh = NA, nm = NA) {
    k <<- k + 1L
    rows[[k]] <<- bind_cols(tibble::tibble(sim = i, condition = cond, generator = gen, rep = rep, n_hc = nh, n_mdd = nm),
                            rank_metrics(ref, rk), truth_metrics(rk, ref, g_main, g_int, g_noise),
                            tibble::tibble(tp50_main_ref = tm_ref$tp50_main, auroc_main_ref = tm_ref$auroc_main))
  }
  add("Reference", "(none)", 0L, ref)
  real0 <- real[real$class == "0", ]; real1 <- real[real$class == "1", ]
  for (g in GENS) {
    syn <- read_syn(g, i); if (is.null(syn)) { msg("sim %02d %s: no synthetic data (skipped)", i, g); next }
    stopifnot(nrow(syn) == 200, setequal(names(syn), names(real))); syn <- syn[, names(real)]
    stopifnot(!anyNA(syn), all(is.finite(as.matrix(syn[, genes]))))
    n0 <- sum(syn$class == "0"); n1 <- sum(syn$class == "1")
    add("Observed", g, 1L, rf_rank(rbind(real, syn), get_seed(i, paste0("rank_obs_", g))), n0, n1)
    set.seed(get_seed(i, paste0("ctrlA_hc_", g)));  a0 <- real0[sample.int(nrow(real0), n0, replace = TRUE), ]
    set.seed(get_seed(i, paste0("ctrlA_mdd_", g))); a1 <- real1[sample.int(nrow(real1), n1, replace = TRUE), ]
    add("Control A", g, 1L, rf_rank(rbind(real, a0, a1), get_seed(i, paste0("ctrlA_rank_", g))), n0, n1)
  }
  for (b in seq_len(N_CTRL_B)) add("Control B", "(none)", b, rf_rank(real, get_seed(i, paste0("ctrlB_rank_", b))))
  for (c in seq_len(N_CTRL_C)) {
    set.seed(get_seed(i, paste0("ctrlC_idx_", c)));  x <- real[sample.int(nrow(real), 200, replace = TRUE), ]
    set.seed(get_seed(i, paste0("ctrlC_perm_", c))); x$class <- sample(x$class)
    add("Control C", "(none)", c, rf_rank(rbind(real, x), get_seed(i, paste0("ctrlC_rank_", c))), sum(x$class == "0"), sum(x$class == "1"))
  }
  msg("sim %02d done: reference Top-50 contains %d of 45 true main-effect genes", i, tm_ref$tp50_main)
}
res <- bind_rows(rows) |> mutate(loss_tp = tp50_main_ref - tp50_main, loss_auroc = auroc_main_ref - auroc_main)
write.csv(res, file.path(res_dir, "per_ranking_metrics.csv"), row.names = FALSE)
nonref <- res |> filter(condition != "Reference")
set.seed(20261028)
boot_ci <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(c(NA_real_, NA_real_))
  bs <- replicate(N_BOOT, mean(sample(x, replace = TRUE))); as.numeric(quantile(bs, c(.025, .975))) }

## ---------- 1) condition summary (mean of replicates within sim, then mean/SD over sims) ----------
KEY <- c("tp50_main", "true_main_dropout", "noise_entries", "auroc_main", "tp50_int", "loss_tp",
         "mean_abs_shift", "median_abs_shift", "dropout_pct", "rbo_p0.95", "jaccard_top50", "spearman")
per_sim <- res |> group_by(sim, condition, generator) |> summarise(across(all_of(KEY), ~ mean(.x, na.rm = TRUE)), n_rankings = n(), .groups = "drop")
cs <- per_sim |> pivot_longer(all_of(KEY), names_to = "metric", values_to = "value") |>
  group_by(condition, generator, metric) |>
  summarise(n_sims = sum(is.finite(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE), .groups = "drop") |>
  arrange(factor(metric, levels = KEY), factor(condition, levels = c("Reference", "Observed", "Control A", "Control B", "Control C")), generator)
write.csv(cs, file.path(res_dir, "table_condition_summary.csv"), row.names = FALSE)

## ---------- 2) paired differences (same simulated data set) ----------
pw <- function(a_cond, a_gen, b_cond, b_gen, label) {
  A <- per_sim |> filter(condition == a_cond, generator == a_gen); B <- per_sim |> filter(condition == b_cond, generator == b_gen)
  J <- inner_join(A, B, by = "sim", suffix = c(".a", ".b"))
  bind_rows(lapply(c("tp50_main", "loss_tp", "true_main_dropout", "noise_entries", "auroc_main", "mean_abs_shift", "dropout_pct", "spearman", "jaccard_top50"), function(m) {
    d <- J[[paste0(m, ".a")]] - J[[paste0(m, ".b")]]; ci <- boot_ci(d)
    tibble::tibble(comparison = label, metric = m, n_sims = sum(is.finite(d)), mean_diff = mean(d, na.rm = TRUE), sd_diff = sd(d, na.rm = TRUE),
                   boot95_low = ci[1], boot95_high = ci[2], n_sims_diff_positive = sum(d > 0, na.rm = TRUE))
  }))
}
pwd <- bind_rows(
  lapply(GENS, function(g) pw("Observed", g, "Control B", "(none)", paste0("Observed ", g, " minus Control B"))),
  lapply(GENS, function(g) pw("Observed", g, "Control A", g, paste0("Observed ", g, " minus Control A"))),
  list(pw("Control C", "(none)", "Control B", "(none)", "Control C minus Control B"),
       pw("Control A", "ARF", "Control B", "(none)", "Control A (ARF-matched) minus Control B"))
)
write.csv(pwd, file.path(res_dir, "table_paired_vs_controls.csv"), row.names = FALSE)

## ---------- 3) which summary tracks the loss of true signal? (within-simulation Spearman) ----------
DISP_M <- c("mean_abs_shift", "median_abs_shift", "dropout_pct")           # larger = more dilution
AGR_M  <- c("spearman", "kendall", "rbo_p0.95", "rbo_p0.98", "jaccard_top25", "jaccard_top50", "jaccard_top100")   # larger = less dilution
orient <- function(x, m) if (m %in% DISP_M) x else -x
tracking <- function(df, outcome, subset_label) {
  bind_rows(lapply(c(DISP_M, AGR_M), function(m) {
    r <- df |> group_by(sim) |> summarise(rho = suppressWarnings(cor(orient(.data[[m]], m), .data[[outcome]], method = "spearman")), .groups = "drop")
    ci <- boot_ci(r$rho)
    tibble::tibble(subset = subset_label, outcome = outcome, metric = m, n_sims = sum(is.finite(r$rho)),
                   mean_rho = mean(r$rho, na.rm = TRUE), boot95_low = ci[1], boot95_high = ci[2])
  }))
}
trk <- bind_rows(tracking(nonref, "loss_tp", "all rankings (incl. Control C)"), tracking(nonref, "loss_auroc", "all rankings (incl. Control C)"),
                 tracking(nonref |> filter(condition != "Control C"), "loss_tp", "without Control C"),
                 tracking(nonref |> filter(condition != "Control C"), "loss_auroc", "without Control C"))
write.csv(trk, file.path(res_dir, "table_metric_vs_truth_loss.csv"), row.names = FALSE)

## ---------- 4) incremental R^2 of dilution summaries beyond agreement / overlap measures ----------
INC_DISP <- c("mean_abs_shift", "median_abs_shift", "dropout_pct"); INC_AGR <- c("spearman", "rbo_p0.95", "jaccard_top50")
r2 <- function(df, outcome, preds) {
  if (!length(preds)) return(0)
  d <- df |> group_by(sim) |> mutate(across(all_of(c(outcome, preds)), ~ .x - mean(.x))) |> ungroup()   # sim fixed effects
  X <- as.matrix(d[, preds]); y <- d[[outcome]]
  if (sd(y) == 0) return(NA_real_); 1 - sum(residuals(lm(y ~ X - 1))^2) / sum(y^2)
}
inc <- function(df, outcome, label) {
  one <- function(d) c(agree = r2(d, outcome, INC_AGR), disp = r2(d, outcome, INC_DISP), both = r2(d, outcome, c(INC_AGR, INC_DISP)))
  est <- one(df); sims_u <- unique(df$sim)
  bs <- t(replicate(N_BOOT, { s <- sample(sims_u, replace = TRUE)
    one(bind_rows(lapply(seq_along(s), function(j) mutate(df[df$sim == s[j], ], sim = j)))) }))
  dd <- cbind(add_disp_to_agree = bs[, "both"] - bs[, "agree"], add_agree_to_disp = bs[, "both"] - bs[, "disp"])
  tibble::tibble(subset = label, outcome = outcome,
    R2_agreement_only = est[["agree"]], R2_dilution_only = est[["disp"]], R2_both = est[["both"]],
    delta_R2_adding_dilution = est[["both"]] - est[["agree"]],
    delta_low = quantile(dd[, 1], .025, na.rm = TRUE), delta_high = quantile(dd[, 1], .975, na.rm = TRUE),
    delta_R2_adding_agreement = est[["both"]] - est[["disp"]],
    delta2_low = quantile(dd[, 2], .025, na.rm = TRUE), delta2_high = quantile(dd[, 2], .975, na.rm = TRUE))
}
msg("Incremental R^2 with %d bootstrap resamples over simulations ...", N_BOOT)
incr <- bind_rows(inc(nonref, "loss_tp", "all rankings (incl. Control C)"), inc(nonref |> filter(condition != "Control C"), "loss_tp", "without Control C"),
                  inc(nonref, "loss_auroc", "all rankings (incl. Control C)"), inc(nonref |> filter(condition != "Control C"), "loss_auroc", "without Control C"))
write.csv(incr, file.path(res_dir, "table_incremental_R2.csv"), row.names = FALSE)

## ---------- 5) can each summary tell a perturbation from ranking noise (Control B)? ----------
det_auc <- function(x, y) mean(outer(x, y, ">")) + 0.5 * mean(outer(x, y, "=="))
det <- bind_rows(lapply(c(DISP_M, AGR_M), function(m) {
  bind_rows(lapply(c("Control C", "Control A", "Observed"), function(cd) {
    a <- sapply(sims, function(i) { x <- nonref |> filter(sim == i, condition == cd) |> pull(.data[[m]]); y <- nonref |> filter(sim == i, condition == "Control B") |> pull(.data[[m]])
      if (!length(x) || !length(y)) NA_real_ else det_auc(orient(x, m), orient(y, m)) })
    tibble::tibble(metric = m, perturbation = cd, n_sims = sum(is.finite(a)), mean_AUC_vs_ControlB = mean(a, na.rm = TRUE))
  }))
}))
write.csv(det, file.path(res_dir, "table_detection_auc.csv"), row.names = FALSE)

## ---------- 6) figures ----------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  nonref$group <- ifelse(nonref$condition == "Observed", nonref$generator, nonref$condition)
  nonref$group <- factor(nonref$group, levels = c("ARF", "CTGAN", "TVAE", "Control A", "Control B", "Control C"))
  ref_tp <- res |> filter(condition == "Reference") |> select(sim, tp50_main)
  p1 <- ggplot(nonref, aes(group, tp50_main)) + geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.4) +
    geom_jitter(aes(colour = group), width = 0.15, height = 0, size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = mean(ref_tp$tp50_main), linetype = 2, colour = "grey40") +
    labs(x = NULL, y = "True main-effect genes in the Top-50", caption = "Dashed line: mean of the reference rankings. Points: one ranking.") +
    theme_bw(base_size = 9) + theme(legend.position = "none", panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(res_dir, "fig_truth_recovery.png"), p1, width = 5, height = 3.5, dpi = 300)
  p2 <- ggplot(nonref, aes(mean_abs_shift, loss_tp, colour = group)) + geom_point(size = 1.2, alpha = 0.7) +
    labs(x = "Mean absolute rank shift of the reference Top-50", y = "Loss of true main-effect genes from the Top-50", colour = NULL) +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(res_dir, "fig_loss_vs_dilution.png"), p2, width = 5.5, height = 3.5, dpi = 300)
}
writeLines(c(sprintf("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf("quick: %s", QUICK),
             sprintf("simulations: %d; RF trees: %d; control B per sim: %d; control C per sim: %d", length(sims), RANK_TREES, N_CTRL_B, N_CTRL_C),
             sprintf("R %s; ranger %s", getRversion(), as.character(packageVersion("ranger")))), file.path(res_dir, "run_info.txt"))

## ---------- 7) console summary ----------
msg("\n=== Condition summary: mean over simulations ===")
print(as.data.frame(cs |> filter(metric %in% c("tp50_main", "loss_tp", "mean_abs_shift", "dropout_pct", "spearman", "jaccard_top50")) |>
      mutate(across(c(mean, sd), ~ round(.x, 3))) |> select(condition, generator, metric, n_sims, mean, sd)), row.names = FALSE)
msg("\n=== Within-simulation Spearman between each summary and the loss of true main-effect genes ===")
print(as.data.frame(trk |> filter(outcome == "loss_tp") |> mutate(across(c(mean_rho, boot95_low, boot95_high), ~ round(.x, 2))) |>
      select(subset, metric, n_sims, mean_rho, boot95_low, boot95_high)), row.names = FALSE)
msg("\n=== Incremental R^2 ===")
print(as.data.frame(incr |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
msg("\nResults written to: %s", res_dir)
