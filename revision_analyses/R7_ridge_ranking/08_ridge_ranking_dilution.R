############################################################
# 08_ridge_ranking_dilution.R   (non-tree sensitivity analysis)
#
# Reviewer 4 (comment 2): "rerun the key analyses with a non-tree evaluator".
# The rank-agreement and dilution analyses of the manuscript rank genes with
# random-forest permutation importance. Here the SAME synthetic draws and the
# SAME controls are re-ranked with a NON-TREE model: ridge logistic regression,
# genes ranked by the absolute standardized coefficient.
#
# Rankings
#   ridge_cv     (main)  predictors standardized with the data being fitted;
#                        ridge penalty chosen by stratified 5-fold CV (lambda.min)
#                        in EVERY fit; CV folds drawn with the seed of that fit.
#   ridge_fixed  (sens.) same, but the penalty is fixed at the value chosen for
#                        the reference ranking, so no ranking differs by lambda.
#   Reference ranking = ridge on the 157 real subjects (CV seed 42).
#
# Conditions (identical to the random-forest analysis 05c_evaluate_repeated_draws.R)
#   Observed   real + synthetic draw d        (ARF, CTGAN, CTGAN_long, TVAE; draws 1..30)
#   Control A  real + 200 rows resampled from the real data, same HC/MDD counts
#   Control B  real data re-ranked with other seeds (here: other CV folds only;
#              ridge has no other source of randomness)
#
# Metrics: Spearman, Kendall, RBO (p = 0.95, 0.98), Top-K Jaccard (25, 50, 100),
# mean / median absolute rank shift and dropout (%) of the reference Top-50.
# Everything is descriptive and conditional on the same 157 real subjects.
#
# Needs (in the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory)): gene_with_label.csv, synthetic_rna_arf.csv /
#   synthetic_ctgan.csv / synthetic_rna_tvae.csv (optional "draw 0"), and the
#   folder reviewer_revision/R4_repeated_draws with draw_seeds.csv and draws/<GEN>/draw_XX.csv
# Writes to the data directory (environment variable SYNTH_EVAL_BASE_DIR, default: current directory)/reviewer_revision/R7_ridge_ranking/
#
# Run (Terminal):   Rscript 08_ridge_ranking_dilution.R
# Smoke test:       RD_QUICK=1 Rscript 08_ridge_ranking_dilution.R   (3 draws, 5 controls; results go to *_QUICK, not usable)
# Time: roughly 20-40 minutes (one cross-validated ridge fit takes a few seconds).
############################################################
base_dir  <- Sys.getenv("RD_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
rd_dir    <- file.path(base_dir, "reviewer_revision", "R4_repeated_draws")
QUICK     <- identical(Sys.getenv("RD_QUICK", "0"), "1")
SFX       <- if (QUICK) "_QUICK" else ""
draws_dir <- file.path(rd_dir, "draws")
out_dir   <- file.path(base_dir, "reviewer_revision", paste0("R7_ridge_ranking", SFX))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MAX_DRAW  <- if (QUICK) 3L else as.integer(Sys.getenv("RD_MAX_DRAW", "30"))
N_CTRL_B  <- if (QUICK) 5L else as.integer(Sys.getenv("RD_CTRL_B", "200"))
REF_SEED  <- 42L; TOP_K <- 50L; N_BOOT <- 2000L; N_FOLDS <- 5L
GENERATORS <- c("ARF", "CTGAN", "CTGAN_long", "TVAE")

suppressPackageStartupMessages({ library(glmnet); library(dplyr); library(tidyr) })
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
t0 <- Sys.time()

## ---------- 1) data and seeds ----------
real <- read.csv(file.path(base_dir, "gene_with_label.csv"), check.names = FALSE)
genes <- setdiff(names(real), "class")
stopifnot(nrow(real) == 157, length(genes) == 1783, setequal(unique(real$class), c(0, 1)),
          all(is.finite(as.matrix(real[, genes]))))
to_class <- function(x) { z <- as.character(x); stopifnot(all(z %in% c("0", "1")))
  factor(ifelse(z == "1", "MDD", "HC"), levels = c("HC", "MDD")) }
real$class <- to_class(real$class)
real_hc <- real[real$class == "HC", ]; real_mdd <- real[real$class == "MDD", ]
f_seeds <- file.path(rd_dir, "draw_seeds.csv")
if (!file.exists(f_seeds)) stop("draw_seeds.csv not found in ", rd_dir, " (run 05_prepare_repeated_draws.R first)")
seeds <- read.csv(f_seeds, stringsAsFactors = FALSE)
get_seed <- function(g, d, role) {
  s <- seeds$seed[seeds$generator == g & seeds$draw == d & seeds$role == role]
  if (length(s) != 1) stop("Seed not found: ", g, " ", d, " ", role); s
}

## ---------- 2) ridge ranking ----------
std_x <- function(df) {                       # standardize with the data being fitted
  x <- as.matrix(df[, genes]); mu <- colMeans(x); s <- apply(x, 2, sd)
  xs <- sweep(sweep(x, 2, mu), 2, ifelse(s > 0, s, 1), "/"); xs[, s == 0] <- 0; xs
}
strat_folds <- function(cls, seed) {          # stratified fold ids drawn with the given seed
  set.seed(seed); id <- integer(length(cls))
  for (lv in levels(cls)) { ix <- which(cls == lv); id[ix] <- sample(rep(seq_len(N_FOLDS), length.out = length(ix))) }
  id
}
ridge_rank <- function(df, seed, lambda = NULL) {
  xs <- std_x(df); y <- as.integer(df$class == "MDD")
  if (is.null(lambda)) {
    cv <- cv.glmnet(xs, y, family = "binomial", alpha = 0, standardize = FALSE, type.measure = "deviance",
                    foldid = strat_folds(df$class, seed))
    lam <- cv$lambda.min; b <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  } else {
    fit <- glmnet(xs, y, family = "binomial", alpha = 0, standardize = FALSE, lambda = lambda)
    lam <- lambda; b <- as.numeric(coef(fit, s = lambda))[-1]
  }
  r <- tibble::tibble(gene = genes, importance = abs(b)) |> arrange(desc(importance), gene) |> mutate(rank = row_number())
  attr(r, "lambda") <- lam; r
}

msg("Reference ridge ranking (real data, CV seed %d) ...", REF_SEED)
ref_cv <- ridge_rank(real, REF_SEED)
LAMBDA_REF <- attr(ref_cv, "lambda")
ref_fx <- ridge_rank(real, REF_SEED, lambda = LAMBDA_REF)
msg("  lambda.min = %.4f (fixed-lambda variant uses the same value)", LAMBDA_REF)
write.csv(ref_cv, file.path(out_dir, "reference_ranking_ridge.csv"), row.names = FALSE)

## ---------- 3) metrics (identical to 05c) ----------
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
    mean_abs_shift = mean(abs(sh)), median_abs_shift = median(abs(sh)), dropout_pct = 100 * mean(rb[top] > K))
}
METRICS <- c("spearman", "kendall", "rbo_p0.95", "rbo_p0.98", "jaccard_top25", "jaccard_top50",
             "jaccard_top100", "mean_abs_shift", "median_abs_shift", "dropout_pct")
DISPLACEMENT <- c("mean_abs_shift", "median_abs_shift", "dropout_pct")

read_draw <- function(f) {
  syn <- read.csv(f, check.names = FALSE)
  if (!setequal(names(syn), names(real))) stop("Columns differ from the real data in ", f)
  syn <- syn[, names(real)]
  if (anyNA(syn) || any(!is.finite(as.matrix(syn[, genes])))) stop("Missing/non-finite values in ", f)
  syn$class <- to_class(syn$class); syn
}
rows <- list(); k <- 0L
add <- function(df) { k <<- k + 1L; rows[[k]] <<- df }
score <- function(ranking, cond, g, d, nh, nm, manu = FALSE) {      # one ranking -> rows for both variants
  for (v in c("ridge_cv", "ridge_fixed")) {
    rk <- if (v == "ridge_cv") ranking$cv else ranking$fx
    ref <- if (v == "ridge_cv") ref_cv else ref_fx
    add(bind_cols(tibble::tibble(variant = v, condition = cond, generator = g, draw = d, manuscript_draw = manu,
                                 n_hc = nh, n_mdd = nm, lambda = attr(rk, "lambda")), rank_metrics(ref, rk)))
  }
}
both <- function(df, seed) list(cv = ridge_rank(df, seed), fx = ridge_rank(df, seed, lambda = LAMBDA_REF))

## ---------- 4) manuscript draws (draw 0), optional ----------
msch <- list(ARF = c("synthetic_rna_arf.csv", "synthetic_arf.csv"),
             CTGAN = c("synthetic_rna_ctgan.csv", "synthetic_ctgan.csv", "synthetic_cgan.csv"),
             TVAE = c("synthetic_rna_tvae.csv", "synthetic_tvae.csv"))
for (g in names(msch)) {
  f <- file.path(base_dir, msch[[g]]); f <- f[file.exists(f)][1]
  if (is.na(f)) { msg("manuscript draw for %s: file not found (skipped)", g); next }
  syn <- read_draw(f)
  score(both(rbind(real, syn[, names(real)]), REF_SEED), "Observed", g, 0L, sum(syn$class == "HC"), sum(syn$class == "MDD"), manu = TRUE)
  msg("manuscript draw %s done", g)
}

## ---------- 5) Control B ----------
msg("Control B: %d re-rankings of the real data (other CV folds) ...", N_CTRL_B)
for (b in seq_len(N_CTRL_B)) {
  score(both(real, get_seed("(none)", b, "ctrlB_rank")), "Control B", "(none)", b, NA_integer_, NA_integer_)
  if (b %% 50 == 0) msg("  control B %d/%d  (%.1f min)", b, N_CTRL_B, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

## ---------- 6) Observed + Control A for each draw ----------
for (g in GENERATORS) {
  files <- file.path(draws_dir, g, sprintf("draw_%02d.csv", seq_len(MAX_DRAW)))
  files <- files[file.exists(files)]
  if (!length(files)) { msg("%s: no draws found", g); next }
  msg("%s: evaluating %d draws ...", g, length(files))
  for (f in files) {
    d <- as.integer(sub("^draw_(\\d+)\\.csv$", "\\1", basename(f)))
    syn <- read_draw(f); nh <- sum(syn$class == "HC"); nm <- sum(syn$class == "MDD"); stopifnot(nh + nm == 200)
    score(both(rbind(real, syn), get_seed(g, d, "rank")), "Observed", g, d, nh, nm)
    set.seed(get_seed(g, d, "ctrlA_hc"));  i1 <- sample.int(nrow(real_hc),  nh, replace = TRUE)
    set.seed(get_seed(g, d, "ctrlA_mdd")); i2 <- sample.int(nrow(real_mdd), nm, replace = TRUE)
    score(both(rbind(real, real_hc[i1, ], real_mdd[i2, ]), get_seed(g, d, "ctrlA_rank")), "Control A", g, d, nh, nm)
    if (d %% 5 == 0) msg("  %s draw %d done (%.1f min)", g, d, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
}
all_res <- bind_rows(rows)
write.csv(all_res, file.path(out_dir, "per_ranking_all_metrics_ridge.csv"), row.names = FALSE)

## ---------- 7) summaries ----------
q <- function(x, p) as.numeric(quantile(x, p, type = 7))
rep_res <- all_res |> filter(!manuscript_draw)
long <- rep_res |> pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value")
summ <- long |> group_by(variant, condition, generator, metric) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), min = min(value), p2.5 = q(value, .025),
            median = median(value), p97.5 = q(value, .975), max = max(value), .groups = "drop") |>
  arrange(variant, factor(metric, levels = METRICS), condition, generator)
write.csv(summ, file.path(out_dir, "summary_by_condition_ridge.csv"), row.names = FALSE)

extreme_count <- function(obs, ctrl, metric) if (metric %in% DISPLACEMENT) sum(ctrl >= obs) else sum(ctrl <= obs)
beyond <- function(obs, ctrl, metric) if (metric %in% DISPLACEMENT) obs > max(ctrl) else obs < min(ctrl)
ex <- list()
for (v in unique(rep_res$variant)) for (g in unique(rep_res$generator[rep_res$condition == "Observed"])) for (mt in METRICS) {
  sel <- rep_res$variant == v
  ob <- rep_res[[mt]][sel & rep_res$condition == "Observed" & rep_res$generator == g]
  cB <- rep_res[[mt]][sel & rep_res$condition == "Control B"]
  cA <- rep_res[[mt]][sel & rep_res$condition == "Control A" & rep_res$generator == g]
  pB <- vapply(ob, function(o) (1 + extreme_count(o, cB, mt)) / (length(cB) + 1), numeric(1))
  ex[[length(ex) + 1]] <- tibble::tibble(
    variant = v, generator = g, metric = mt, n_draws = length(ob), observed_mean = mean(ob),
    controlA_mean = mean(cA), controlB_mean = mean(cB),
    n_draws_beyond_all_controlB = sum(vapply(ob, beyond, logical(1), ctrl = cB, metric = mt)),
    n_draws_beyond_all_controlA = sum(vapply(ob, beyond, logical(1), ctrl = cA, metric = mt)),
    emp_exceedance_controlB_median = median(pB), emp_exceedance_controlB_min_possible = 1 / (length(cB) + 1))
}
write.csv(bind_rows(ex), file.path(out_dir, "exceedance_vs_controls_ridge.csv"), row.names = FALSE)

set.seed(20261029)
gens_present <- intersect(GENERATORS, unique(rep_res$generator[rep_res$condition == "Observed"]))
cmp <- list()
if (length(gens_present) >= 2) for (v in unique(rep_res$variant)) for (pr in combn(gens_present, 2, simplify = FALSE)) for (mt in METRICS) {
  sel <- rep_res$variant == v & rep_res$condition == "Observed"
  x <- rep_res[[mt]][sel & rep_res$generator == pr[1]]; y <- rep_res[[mt]][sel & rep_res$generator == pr[2]]
  if (length(x) < 3 || length(y) < 3) next
  bs <- replicate(N_BOOT, mean(sample(x, replace = TRUE)) - mean(sample(y, replace = TRUE)))
  cmp[[length(cmp) + 1]] <- tibble::tibble(variant = v, first = pr[1], second = pr[2], metric = mt, n_first = length(x), n_second = length(y),
    mean_first = mean(x), mean_second = mean(y), mean_diff_first_minus_second = mean(x) - mean(y), boot95_low = q(bs, .025), boot95_high = q(bs, .975))
}
write.csv(bind_rows(cmp), file.path(out_dir, "generator_pairwise_comparison_ridge.csv"), row.names = FALSE)

fmt <- function(x, dg) sprintf("%.*f (%.*f) [%.*f to %.*f]", dg, mean(x), dg, sd(x), dg, min(x), dg, max(x))
tab_metrics <- c(mean_abs_shift = 1, median_abs_shift = 1, dropout_pct = 1, `rbo_p0.95` = 3, jaccard_top50 = 3, spearman = 3)
tb <- rep_res |> group_by(variant, condition, generator) |> group_modify(function(df, key) {
  out <- tibble::tibble(n = nrow(df)); for (mt in names(tab_metrics)) out[[mt]] <- fmt(df[[mt]], tab_metrics[[mt]]); out }) |> ungroup() |>
  arrange(variant, factor(condition, levels = c("Observed", "Control A", "Control B")), generator)
write.csv(tb, file.path(out_dir, "table_mean_sd_range_ridge.csv"), row.names = FALSE)

# how does the ridge reference ranking relate to the random-forest reference ranking (if available)?
f_rf <- file.path(rd_dir, "results", "reference_ranking.csv")
rf_note <- "random-forest reference ranking not found (skipped)"
if (file.exists(f_rf)) {
  rf <- read.csv(f_rf); rf$rank <- rank(-rf$importance, ties.method = "first")
  rf <- rf[match(ref_cv$gene, rf$gene), ]
  top_r <- ref_cv$gene[ref_cv$rank <= TOP_K]; top_f <- rf$gene[rf$rank <= TOP_K]
  rf_note <- sprintf("ridge vs random-forest reference ranking: Spearman = %.3f; Top-50 overlap = %d genes (Jaccard %.3f)",
                     cor(ref_cv$rank, rf$rank, method = "spearman"), length(intersect(top_r, top_f)),
                     length(intersect(top_r, top_f)) / length(union(top_r, top_f)))
}
# in-sample fit of the reference ridge (for context only; performance on held-out subjects is in Table 1)
xs <- std_x(real); fit <- glmnet(xs, as.integer(real$class == "MDD"), family = "binomial", alpha = 0, standardize = FALSE, lambda = LAMBDA_REF)
writeLines(c(
  sprintf("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf("quick: %s", QUICK),
  sprintf("reference ridge: lambda.min = %.5f (CV seed %d, %d stratified folds); Top-K = %d; control B replicates = %d", LAMBDA_REF, REF_SEED, N_FOLDS, TOP_K, N_CTRL_B),
  rf_note,
  sprintf("draws per generator (observed): %s", paste(names(table(rep_res$generator[rep_res$condition == "Observed" & rep_res$variant == "ridge_cv"])),
          table(rep_res$generator[rep_res$condition == "Observed" & rep_res$variant == "ridge_cv"]), sep = "=", collapse = ", ")),
  sprintf("seed file md5: %s", unname(tools::md5sum(f_seeds))),
  sprintf("R %s; glmnet %s", getRversion(), as.character(packageVersion("glmnet")))), file.path(out_dir, "run_info_ridge.txt"))

## ---------- 8) console summary ----------
msg("\n%s", rf_note)
msg("\n================ SUMMARY (ridge_cv): mean (SD) [min to max] ================")
print(as.data.frame(tb |> filter(variant == "ridge_cv") |> select(condition, generator, n, mean_abs_shift, dropout_pct, `rbo_p0.95`, jaccard_top50, spearman)), right = FALSE)
msg("\n---- Observed draws beyond ALL Control B / Control A values (mean absolute rank shift), ridge_cv ----")
print(as.data.frame(bind_rows(ex) |> filter(variant == "ridge_cv", metric == "mean_abs_shift") |>
        select(generator, n_draws, observed_mean, controlA_mean, controlB_mean, n_draws_beyond_all_controlA, n_draws_beyond_all_controlB)), row.names = FALSE)
msg("\nResults written to: %s  (%.1f min)", out_dir, as.numeric(difftime(Sys.time(), t0, units = "mins")))
