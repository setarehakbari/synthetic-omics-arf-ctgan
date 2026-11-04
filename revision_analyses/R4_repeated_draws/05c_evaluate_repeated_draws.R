############################################################
# 05c_evaluate_repeated_draws.R   (repeated-draws pipeline, step 3)
#
# Evaluates the repeated synthetic-data draws (steps 1-2) with the SAME
# ranking pipeline and metrics as the main analysis (Table 2):
#   reference ranking = random forest on the real data, permutation
#   importance, 3000 trees, seed 42; every augmented ranking = the same
#   random forest on real + 200 synthetic subjects.
#
# For every generator and draw:
#   Observed   real + synthetic draw d
#   Control A  real + 200 rows resampled (with replacement) from the real
#              data, same HC/MDD counts as that synthetic draw (no generator)
# Once (not tied to a generator):
#   Control B  real data only, re-ranked with N_CTRL_B other RF seeds
#              (pure random-forest noise)
#
# Metrics (all against the fixed reference ranking)
#   Spearman, Kendall, RBO (p = 0.95 as in Table 2; p = 0.98 as sensitivity),
#   Top-K Jaccard (K = 25, 50, 100), mean / median absolute rank shift and
#   % dropout of the real Top-50, mean / median signed shift.
#
# Reported: mean, SD, min, 2.5/50/97.5 percentiles, max over draws; how
# many draws lie beyond every control value; empirical p-values
# (1 + #controls as or more extreme) / (B + 1); and descriptive
# pairwise generator comparisons (mean difference with bootstrap 95% CI,
# probability of superiority, Wilcoxon rank-sum p-value).
#
# All draws are conditional on the same 157 real subjects: the intervals
# describe how much a result varies from one run of a generator to the next
# and how it compares with ranking noise, NOT sampling variability of the
# real population. Everything is descriptive.
#
# If the original single draws of the manuscript are found
# (synthetic_rna_arf.csv, synthetic_rna_ctgan.csv, synthetic_rna_tvae.csv
# in the base folder) they are re-evaluated as "draw 0" (seed 42, as in
# Table 2) to check that this script reproduces Table 2 and to show where
# the manuscript draw lies within the repeated draws.
#
# Usage
#   Rscript 05c_evaluate_repeated_draws.R
#   RD_MAX_DRAW=30 (default) uses draws 1..30 of each generator that exist.
#   RD_QUICK=1 Rscript 05c_evaluate_repeated_draws.R     # smoke test only
############################################################
base_dir <- Sys.getenv("RD_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
rd_dir   <- file.path(base_dir, "reviewer_revision", "R4_repeated_draws")
QUICK    <- identical(Sys.getenv("RD_QUICK", "0"), "1")
SFX      <- if (QUICK) "_QUICK" else ""
draws_dir <- file.path(rd_dir, paste0("draws", SFX))
res_dir   <- file.path(rd_dir, paste0("results", SFX))
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
MAX_DRAW   <- as.integer(Sys.getenv("RD_MAX_DRAW", "30"))
N_CTRL_B   <- if (QUICK) 5L else as.integer(Sys.getenv("RD_CTRL_B", "200"))
RANK_TREES <- if (QUICK) 300L else 3000L
REF_SEED   <- 42L
TOP_K      <- 50L
N_BOOT     <- 2000L
GENERATORS <- c("ARF", "CTGAN", "CTGAN_long", "TVAE")

suppressPackageStartupMessages({ library(ranger); library(dplyr); library(tidyr) })
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ---------- 1) Real data, seeds, reference ranking ----------
real <- read.csv(file.path(base_dir, "gene_with_label.csv"), check.names = FALSE)
genes <- setdiff(names(real), "class")
stopifnot(nrow(real) == 157, length(genes) == 1783, setequal(unique(real$class), c(0, 1)),
          all(is.finite(as.matrix(real[, genes]))))
to_class <- function(x) {
  z <- as.character(x); stopifnot(all(z %in% c("0", "1")))
  factor(ifelse(z == "1", "MDD", "HC"), levels = c("HC", "MDD"))
}
real$class <- to_class(real$class)
real_hc  <- real[real$class == "HC", ];  real_mdd <- real[real$class == "MDD", ]
seeds <- read.csv(file.path(rd_dir, "draw_seeds.csv"), stringsAsFactors = FALSE)
get_seed <- function(g, d, role) {
  s <- seeds$seed[seeds$generator == g & seeds$draw == d & seeds$role == role]
  if (length(s) != 1) stop("Seed not found: ", g, " ", d, " ", role); s
}

rf_rank <- function(df, seed) {
  m <- ranger(class ~ ., data = df, num.trees = RANK_TREES, importance = "permutation",
              scale.permutation.importance = TRUE, seed = seed)
  tibble::tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |> mutate(rank = row_number())
}

msg("Reference ranking (real data only, seed %d, %d trees) ...", REF_SEED, RANK_TREES)
ref <- rf_rank(real, REF_SEED)
write.csv(ref, file.path(res_dir, "reference_ranking.csv"), row.names = FALSE)

## ---------- 2) Metrics ----------
# RBO with extrapolation, identical to the definition in Table 2:
#   (1-p) * sum_{d=1..k} (X_d/d) p^(d-1) + (X_k/k) p^k,  X_d = overlap of the two top-d lists.
# A gene is in both top-d lists exactly when max(rank_A, rank_B) <= d, hence X_d is a cumulative count.
rbo_fast <- function(ra, rb, p) {
  k <- length(ra); X <- cumsum(tabulate(pmax(ra, rb), nbins = k)); d <- seq_len(k); agree <- X / d
  (1 - p) * sum(agree * p^(d - 1)) + agree[k] * p^k
}
rank_metrics <- function(ref, aug, K = TOP_K) {
  ra <- ref$rank; rb <- aug$rank[match(ref$gene, aug$gene)]
  stopifnot(!anyNA(rb), length(rb) == length(ra))
  jac <- function(k) { a <- ref$gene[ra <= k]; b <- ref$gene[rb <= k]; length(intersect(a, b)) / length(union(a, b)) }
  top <- ra <= K; sh <- rb[top] - ra[top]
  tibble::tibble(
    spearman = cor(ra, rb, method = "spearman"), kendall = cor(ra, rb, method = "kendall"),
    `rbo_p0.95` = rbo_fast(ra, rb, 0.95), `rbo_p0.98` = rbo_fast(ra, rb, 0.98),
    jaccard_top25 = jac(25), jaccard_top50 = jac(50), jaccard_top100 = jac(100),
    mean_abs_shift = mean(abs(sh)), median_abs_shift = median(abs(sh)), dropout_pct = 100 * mean(rb[top] > K),
    mean_signed_shift = mean(sh), median_signed_shift = median(sh))
}
METRICS      <- c("spearman", "kendall", "rbo_p0.95", "rbo_p0.98", "jaccard_top25", "jaccard_top50",
                  "jaccard_top100", "mean_abs_shift", "median_abs_shift", "dropout_pct")
DISPLACEMENT <- c("mean_abs_shift", "median_abs_shift", "dropout_pct")      # larger = more dilution
SIGNED       <- c("mean_signed_shift", "median_signed_shift")
# self-test of the fast RBO against the plain definition used in Table2_rank_metrics_with_TVAE.R
rbo_plain <- function(l1, l2, p) { k <- min(length(l1), length(l2))
  ov <- vapply(seq_len(k), function(d) length(intersect(l1[1:d], l2[1:d])) / d, numeric(1))
  (1 - p) * sum(ov * p^(0:(k - 1))) + ov[k] * p^k }
rbo_check <- { a <- ref$gene; b <- sample(a); ra <- seq_along(a); rb <- match(a, b)
  abs(rbo_plain(a, b, 0.95) - rbo_fast(ra, rb, 0.95)) < 1e-10 }
if (!rbo_check) stop("fast RBO does not match the plain definition")
msg("RBO implementation self-test passed.")

read_draw <- function(f) {
  syn <- read.csv(f, check.names = FALSE)
  if (!setequal(names(syn), names(real))) stop("Columns differ from the real data in ", f)
  syn <- syn[, names(real)]
  if (anyNA(syn) || any(!is.finite(as.matrix(syn[, genes])))) stop("Missing/non-finite values in ", f)
  syn$class <- to_class(syn$class); syn
}

## ---------- 3) Manuscript draws (draw 0) ----------
msch <- list(
  ARF   = c("synthetic_rna_arf.csv", "synthetic_arf.csv"),
  CTGAN = c("synthetic_rna_ctgan.csv", "synthetic_ctgan.csv", "synthetic_cgan.csv"),
  TVAE  = c("synthetic_rna_tvae.csv", "synthetic_tvae.csv"))
rows <- list(); k <- 0L
add <- function(df) { k <<- k + 1L; rows[[k]] <<- df }
expected_T2 <- c(ARF = 418.36, CTGAN = 481.66, TVAE = 413.38)    # Table 2 of the manuscript (mean abs. rank shift)
for (g in names(msch)) {
  f <- file.path(base_dir, msch[[g]]); f <- f[file.exists(f)][1]
  if (is.na(f)) { msg("manuscript draw for %s: file not found (skipped)", g); next }
  syn <- read_draw(f)
  m <- rank_metrics(ref, rf_rank(rbind(real, syn[, names(real)]), REF_SEED))
  add(bind_cols(tibble::tibble(condition = "Observed", generator = g, draw = 0L, manuscript_draw = TRUE,
                               n_hc = sum(syn$class == "HC"), n_mdd = sum(syn$class == "MDD")), m))
  if (!QUICK) msg("manuscript draw %s: mean abs. shift = %.2f (Table 2: %.2f) -> %s", g, m$mean_abs_shift,
                  expected_T2[[g]], if (abs(m$mean_abs_shift - expected_T2[[g]]) < 0.01) "REPRODUCED" else "differs from Table 2 (ranger version or file version? please report)")
}

## ---------- 4) Control B (real data, other RF seeds) ----------
msg("Control B: %d re-rankings of the real data ...", N_CTRL_B)
for (b in seq_len(N_CTRL_B)) {
  m <- rank_metrics(ref, rf_rank(real, get_seed("(none)", b, "ctrlB_rank")))
  add(bind_cols(tibble::tibble(condition = "Control B", generator = "(none)", draw = b, manuscript_draw = FALSE,
                               n_hc = NA_integer_, n_mdd = NA_integer_), m))
  if (b %% 50 == 0) msg("  control B %d/%d", b, N_CTRL_B)
}

## ---------- 5) Observed + Control A for every available draw ----------
for (g in GENERATORS) {
  files <- file.path(draws_dir, g, sprintf("draw_%02d.csv", seq_len(MAX_DRAW)))
  files <- files[file.exists(files)]
  if (!length(files)) { msg("%s: no draws found", g); next }
  msg("%s: evaluating %d draws ...", g, length(files))
  for (f in files) {
    d <- as.integer(sub("^draw_(\\d+)\\.csv$", "\\1", basename(f)))
    syn <- read_draw(f); nh <- sum(syn$class == "HC"); nm <- sum(syn$class == "MDD")
    stopifnot(nh + nm == 200)
    m <- rank_metrics(ref, rf_rank(rbind(real, syn), get_seed(g, d, "rank")))
    add(bind_cols(tibble::tibble(condition = "Observed", generator = g, draw = d, manuscript_draw = FALSE,
                                 n_hc = nh, n_mdd = nm), m))
    set.seed(get_seed(g, d, "ctrlA_hc"));  i1 <- sample.int(nrow(real_hc),  nh, replace = TRUE)
    set.seed(get_seed(g, d, "ctrlA_mdd")); i2 <- sample.int(nrow(real_mdd), nm, replace = TRUE)
    mA <- rank_metrics(ref, rf_rank(rbind(real, real_hc[i1, ], real_mdd[i2, ]), get_seed(g, d, "ctrlA_rank")))
    add(bind_cols(tibble::tibble(condition = "Control A", generator = g, draw = d, manuscript_draw = FALSE,
                                 n_hc = nh, n_mdd = nm), mA))
    if (d %% 5 == 0) msg("  %s draw %d done", g, d)
  }
}
all_res <- bind_rows(rows)
write.csv(all_res, file.path(res_dir, "per_draw_all_metrics.csv"), row.names = FALSE)

## ---------- 6) Summaries ----------
rep_res <- all_res |> filter(!manuscript_draw)
long <- rep_res |> pivot_longer(all_of(c(METRICS, SIGNED)), names_to = "metric", values_to = "value")
q <- function(x, p) as.numeric(quantile(x, p, type = 7))
summ <- long |> group_by(condition, generator, metric) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), min = min(value), p2.5 = q(value, .025),
            median = median(value), p97.5 = q(value, .975), max = max(value), .groups = "drop") |>
  arrange(factor(metric, levels = c(METRICS, SIGNED)), condition, generator)
write.csv(summ, file.path(res_dir, "summary_by_condition.csv"), row.names = FALSE)

# class balance of the synthetic draws
cb <- rep_res |> filter(condition == "Observed") |> group_by(generator) |>
  summarise(n_draws = n(), hc_mean = mean(n_hc), hc_min = min(n_hc), hc_max = max(n_hc),
            mdd_mean = mean(n_mdd), mdd_min = min(n_mdd), mdd_max = max(n_mdd), .groups = "drop")
write.csv(cb, file.path(res_dir, "synthetic_class_balance.csv"), row.names = FALSE)

# exceedance / empirical p-values (direction-aware)
extreme_count <- function(obs, ctrl, metric) if (metric %in% DISPLACEMENT) sum(ctrl >= obs) else sum(ctrl <= obs)
beyond <- function(obs, ctrl, metric) if (metric %in% DISPLACEMENT) obs > max(ctrl) else obs < min(ctrl)
beyond975 <- function(obs, ctrl, metric) if (metric %in% DISPLACEMENT) obs > q(ctrl, .975) else obs < q(ctrl, .025)
ex <- list()
for (g in unique(rep_res$generator[rep_res$condition == "Observed"])) for (mt in METRICS) {
  ob <- rep_res[[mt]][rep_res$condition == "Observed"  & rep_res$generator == g]
  cB <- rep_res[[mt]][rep_res$condition == "Control B"]
  cA <- rep_res[[mt]][rep_res$condition == "Control A" & rep_res$generator == g]
  pB <- vapply(ob, function(o) (1 + extreme_count(o, cB, mt)) / (length(cB) + 1), numeric(1))
  pA <- vapply(ob, function(o) (1 + extreme_count(o, cA, mt)) / (length(cA) + 1), numeric(1))
  least <- if (mt %in% DISPLACEMENT) min(ob) else max(ob)               # least-diluted observed draw
  most  <- function(x) if (mt %in% DISPLACEMENT) max(x) else min(x)     # most extreme control value
  ex[[length(ex) + 1]] <- tibble::tibble(
    generator = g, metric = mt, n_draws = length(ob),
    observed_mean = mean(ob), observed_least_diluted_draw = least,
    controlA_most_extreme = most(cA), controlB_most_extreme = most(cB),
    n_draws_beyond_all_controlB = sum(vapply(ob, beyond, logical(1), ctrl = cB, metric = mt)),
    n_draws_beyond_all_controlA = sum(vapply(ob, beyond, logical(1), ctrl = cA, metric = mt)),
    n_draws_beyond_controlB_97.5pct = sum(vapply(ob, beyond975, logical(1), ctrl = cB, metric = mt)),
    emp_p_controlB_median = median(pB), emp_p_controlB_max = max(pB), emp_p_controlB_min_possible = 1 / (length(cB) + 1),
    emp_p_controlA_median = median(pA), emp_p_controlA_max = max(pA), emp_p_controlA_min_possible = 1 / (length(cA) + 1))
}
write.csv(bind_rows(ex), file.path(res_dir, "exceedance_vs_controls.csv"), row.names = FALSE)

# pairwise generator comparison (descriptive; draws are independent given the real data)
set.seed(20261027)
gens_present <- intersect(c("ARF", "CTGAN", "TVAE", "CTGAN_long"), unique(rep_res$generator[rep_res$condition == "Observed"]))
pairs <- if (length(gens_present) >= 2) combn(gens_present, 2, simplify = FALSE) else list()
cmp <- list()
for (pr in pairs) for (mt in METRICS) {
  x <- rep_res[[mt]][rep_res$condition == "Observed" & rep_res$generator == pr[1]]
  y <- rep_res[[mt]][rep_res$condition == "Observed" & rep_res$generator == pr[2]]
  if (length(x) < 5 || length(y) < 5) next
  bs <- replicate(N_BOOT, mean(sample(x, replace = TRUE)) - mean(sample(y, replace = TRUE)))
  wp <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
  cmp[[length(cmp) + 1]] <- tibble::tibble(
    first = pr[1], second = pr[2], metric = mt, n_first = length(x), n_second = length(y),
    mean_first = mean(x), mean_second = mean(y), mean_diff_first_minus_second = mean(x) - mean(y),
    boot95_low = q(bs, .025), boot95_high = q(bs, .975),
    prob_first_greater = mean(outer(x, y, ">")) + 0.5 * mean(outer(x, y, "==")), wilcoxon_p = wp)
}
write.csv(bind_rows(cmp), file.path(res_dir, "generator_pairwise_comparison.csv"), row.names = FALSE)

# compact table for the manuscript / response letter
fmt <- function(x, dg) sprintf("%.*f (%.*f) [%.*f to %.*f]", dg, mean(x), dg, sd(x), dg, min(x), dg, max(x))
tab_metrics <- c(mean_abs_shift = 1, median_abs_shift = 1, dropout_pct = 1, `rbo_p0.95` = 3, jaccard_top50 = 3, spearman = 3)
lab <- function(cond, g) if (cond == "Observed") paste0("Real + ", g) else if (cond == "Control A") paste0("Control A (", g, "-matched)") else "Control B (RF seed only)"
tb <- rep_res |> group_by(condition, generator) |> group_modify(function(df, key) {
  out <- tibble::tibble(label = lab(key$condition, key$generator), n = nrow(df))
  for (mt in names(tab_metrics)) out[[mt]] <- fmt(df[[mt]], tab_metrics[[mt]]); out }) |> ungroup()
tb <- tb |> arrange(factor(condition, levels = c("Observed", "Control A", "Control B")), generator)
ms <- all_res |> filter(manuscript_draw) |> select(generator, all_of(names(tab_metrics)))
if (nrow(ms)) tb$manuscript_single_draw_mean_abs_shift <- ifelse(tb$condition == "Observed",
      round(ms$mean_abs_shift[match(tb$generator, ms$generator)], 1), NA)
write.csv(tb, file.path(res_dir, "table_mean_sd_range.csv"), row.names = FALSE)

# generation logs of all draws (training budgets, run times, versions)
logs <- unlist(lapply(GENERATORS, function(g) list.files(file.path(draws_dir, g), "^log_\\d+\\.csv$", full.names = TRUE)))
if (length(logs)) write.csv(bind_rows(lapply(logs, read.csv, colClasses = "character")), file.path(res_dir, "generation_log_all.csv"), row.names = FALSE)

writeLines(c(
  sprintf("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf("quick: %s", QUICK),
  sprintf("rank trees: %d; reference seed: %d; Top-K: %d; control B replicates: %d", RANK_TREES, REF_SEED, TOP_K, N_CTRL_B),
  sprintf("draws per generator (observed): %s", paste(names(table(rep_res$generator[rep_res$condition == "Observed"])),
                                                     table(rep_res$generator[rep_res$condition == "Observed"]), sep = "=", collapse = ", ")),
  sprintf("seed file md5: %s", unname(tools::md5sum(file.path(rd_dir, "draw_seeds.csv")))),
  sprintf("R %s; ranger %s", getRversion(), as.character(packageVersion("ranger")))),
  file.path(res_dir, "run_info.txt"))

## ---------- 7) Figure ----------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  col <- c(ARF = "#2a78d6", CTGAN = "#eb6834", TVAE = "#1baf7a", CTGAN_long = "#e87ba4")
  lvl <- c(intersect(c("ARF", "CTGAN", "CTGAN_long", "TVAE"), unique(rep_res$generator[rep_res$condition == "Observed"])),
           "Control A", "Control B")
  pd <- rep_res |>
    mutate(group = ifelse(condition == "Observed", generator, condition)) |>
    select(group, condition, generator, all_of(c("mean_abs_shift", "dropout_pct", "rbo_p0.95", "jaccard_top50"))) |>
    pivot_longer(c(mean_abs_shift, dropout_pct, `rbo_p0.95`, jaccard_top50), names_to = "metric", values_to = "value") |>
    mutate(group = factor(group, levels = lvl),
           metric = factor(metric, levels = c("mean_abs_shift", "dropout_pct", "rbo_p0.95", "jaccard_top50"),
                           labels = c("Mean absolute rank shift (Top-50)", "Top-50 dropout (%)", "RBO (p = 0.95)", "Top-50 Jaccard")),
           fill_key = ifelse(condition == "Observed", as.character(generator), "ctrl"))
  mp <- all_res |> filter(manuscript_draw) |>
    pivot_longer(c(mean_abs_shift, dropout_pct, `rbo_p0.95`, jaccard_top50), names_to = "metric", values_to = "value") |>
    filter(generator %in% lvl) |>
    mutate(group = factor(generator, levels = lvl),
           metric = factor(metric, levels = c("mean_abs_shift", "dropout_pct", "rbo_p0.95", "jaccard_top50"),
                           labels = c("Mean absolute rank shift (Top-50)", "Top-50 dropout (%)", "RBO (p = 0.95)", "Top-50 Jaccard")))
  p <- ggplot(pd, aes(group, value)) +
    geom_boxplot(aes(colour = fill_key), width = 0.55, linewidth = 0.4, outlier.shape = NA, fill = NA) +
    geom_jitter(aes(colour = fill_key), width = 0.14, height = 0, size = 1.3, alpha = 0.75) +
    geom_point(data = mp, shape = 23, size = 2.6, stroke = 0.6, fill = "white", colour = "black") +
    scale_colour_manual(values = c(col, ctrl = "#52514e"), guide = "none") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = NULL,
         caption = "Each point is one independent draw (Control A: resampled real rows; Control B: random-forest seed only).\nDiamonds: the single draw reported in the manuscript.") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0, colour = "#52514e"))
  ggsave(file.path(res_dir, "fig_repeated_draws.png"), p, width = 10, height = 3.6, dpi = 300)
  ggsave(file.path(res_dir, "fig_repeated_draws.pdf"), p, width = 10, height = 3.6)
}

## ---------- 8) Console summary ----------
msg("\n================ SUMMARY: mean (SD) [min to max] ================")
print(as.data.frame(tb |> select(label, n, mean_abs_shift, dropout_pct, `rbo_p0.95`, jaccard_top50)), right = FALSE)
if (length(ex)) {
msg("\n---- Observed draws beyond ALL control values (mean absolute rank shift) ----")
print(as.data.frame(bind_rows(ex) |> filter(metric == "mean_abs_shift") |>
        select(generator, n_draws, observed_least_diluted_draw, controlA_most_extreme, controlB_most_extreme,
               n_draws_beyond_all_controlA, n_draws_beyond_all_controlB, emp_p_controlB_max)), row.names = FALSE)
}
msg("\nResults written to: %s", res_dir)
