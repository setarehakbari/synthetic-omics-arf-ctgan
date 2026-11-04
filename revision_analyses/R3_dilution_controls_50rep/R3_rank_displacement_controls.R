############################################################
# R3_rank_displacement_controls.R
#
# Supplementary analysis for the second major-revision round of
# "Multidimensional Evaluation of Synthetic Gene Expression Data with
#  a Dilution Metric for Feature Ranking Preservation" (BMC
#  Bioinformatics).
#
# Addresses Reviewer 3, Major Comment 1:
#   "Gene ranking is not very stable even within the real data itself.
#    How do you know the rank displacement observed after augmenting
#    with synthetic data reflects the synthetic data itself, rather
#    than ordinary sampling variability or correlated-gene rank
#    swapping that would occur anyway?"
#
# This script builds two controls that use ONLY the real data (no
# generator involved), expressed in the SAME units as the manuscript's
# Dilution metric (mean/median absolute rank shift among genes that
# start in the Real Top-K; see Table 2), so they are directly
# comparable to the observed Real-vs-Real+{ARF,CTGAN,TVAE} numbers:
#
#   CONTROL A - Resampling control (no generator):
#     Instead of adding 200 synthetic rows, add 200 rows resampled
#     WITH REPLACEMENT from the real data itself, matching the HC/MDD
#     class balance actually produced by each generator. Refit the RF
#     ranking on Real + resampled-Real and measure displacement
#     relative to the Real-only ranking, repeated B times.
#     NOTE (stated explicitly in the output): resampled rows are
#     repeats of existing individuals, not new independent samples;
#     this control isolates the effect of RE-WEIGHTING the training
#     set, not the effect of genuinely new data.
#
#   CONTROL B - RF stochasticity control (no data change at all):
#     Refit the Real-only RF ranking B times with different random
#     seeds on the SAME 157 rows, and measure displacement against a
#     fixed reference ranking (seed = SEED). This isolates pure
#     algorithmic randomness in permutation importance.
#
#   SUPPLEMENTARY CHECK - Correlated-gene substitution:
#     For each gene that leaves the Real Top-K after augmentation,
#     check whether its most-correlated partner gene (correlation
#     computed from REAL data only, so the check is not circular)
#     entered/stayed in the Top-K in its place. Reported as descriptive
#     evidence, not as an exact decomposition of "signal retained."
#
# USAGE
#   Rscript R3_rank_displacement_controls.R
#
# INPUT (only the first file is required):
#   gene_with_label.csv        - canonical real dataset (157 x 1784)
#   synthetic_rna_arf.csv      - optional; used to read the exact
#   synthetic_ctgan.csv          HC/MDD class balance each generator
#   synthetic_rna_tvae.csv       produced (for Control A), and, if
#                                 present, to also report the ACTUAL
#                                 observed dilution for a side-by-side
#                                 comparison. If a file is missing, its
#                                 generator's row in Control A falls
#                                 back to a documented default class
#                                 balance (the one in the archived
#                                 manuscript files), with a warning.
#
# OUTPUT (written under <base_dir>/reviewer_revision/R3_major1_controls/):
#   controlA_resampling_replicates.csv
#   controlB_seed_stability_replicates.csv
#   correlated_gene_substitution.csv
#   R3_major1_controls_summary.csv   <- the table to quote in the
#                                        response letter
#
# This script never modifies any existing manuscript output file.
############################################################

## =========================================================
## 0) USER SETTINGS
## =========================================================
QUICK_TEST  <- FALSE   # set TRUE for a fast sanity check (fewer trees/reps)
SEED        <- 42
TOP_K       <- 50
RANK_TREES  <- 3000    # matches the RF ranking settings used for Table 2
N_INJECT    <- 200     # matches the number of synthetic rows per generator
B_CONTROL_A <- 50       # resampling-control replicates (per generator)
B_CONTROL_B <- 50       # RF-seed-stochasticity replicates

if (QUICK_TEST) {
  RANK_TREES  <- 300
  B_CONTROL_A <- 3
  B_CONTROL_B <- 3
}

base_dir  <- Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = ".")
real_path <- file.path(base_dir, "gene_with_label.csv")
arf_path  <- file.path(base_dir, "synthetic_rna_arf.csv")
ctg_path  <- file.path(base_dir, "synthetic_ctgan.csv")
tvae_path <- file.path(base_dir, "synthetic_rna_tvae.csv")

out_dir <- file.path(base_dir, "reviewer_revision", "R3_major1_controls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 1) PACKAGES
## =========================================================
required_pkgs <- c("dplyr", "tidyr", "readr", "tibble", "purrr", "ranger")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(", paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
  )
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
  library(ranger)
})

## =========================================================
## 2) HELPERS (kept consistent with evaluate_all_generators.R)
## =========================================================
msg <- function(...) cat(sprintf(...), "\n")

fix_class <- function(x, dataset_name) {
  z <- trimws(as.character(x))
  uz <- unique(z)
  if (all(uz %in% c("0", "1"))) return(factor(z, levels = c("0", "1"), labels = c("HC", "MDD")))
  z_up <- toupper(z)
  if (all(unique(z_up) %in% c("HC", "MDD"))) return(factor(z_up, levels = c("HC", "MDD")))
  stop("Unexpected class labels in ", dataset_name, ": ", paste(sort(unique(z)), collapse = ", "))
}

force_numeric_finite <- function(df, genes, dataset_name) {
  for (g in genes) df[[g]] <- suppressWarnings(as.numeric(df[[g]]))
  if (anyNA(df[, genes, drop = FALSE])) stop("Missing/non-numeric values in ", dataset_name, ".")
  m <- as.matrix(df[, genes, drop = FALSE])
  if (any(!is.finite(m))) stop("Inf/-Inf detected in ", dataset_name, ".")
  df
}

rf_rank <- function(df, seed_value, trees = RANK_TREES) {
  m <- ranger(
    class ~ ., data = df, num.trees = trees,
    importance = "permutation", scale.permutation.importance = TRUE,
    seed = seed_value
  )
  tibble(gene = names(m$variable.importance), importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance), gene) |>
    mutate(rank = row_number())
}

## Dilution metric, in the same units/definition used for Table 2:
## mean/median ABSOLUTE RANK SHIFT among genes that are in the Real
## Top-K, plus the % of those genes that fall out of the Top-K.
dilution_vs_real <- function(rank_real, rank_aug, K = TOP_K) {
  merged <- rank_real |>
    select(gene, rank_real = rank) |>
    inner_join(rank_aug |> select(gene, rank_aug = rank), by = "gene")
  top_real <- rank_real$gene[seq_len(K)]
  detail <- merged |>
    filter(gene %in% top_real) |>
    mutate(
      signed_rank_shift = rank_aug - rank_real,
      absolute_rank_shift = abs(signed_rank_shift),
      dropped = rank_aug > K
    )
  tibble(
    mean_absolute_rank_shift   = mean(detail$absolute_rank_shift),
    median_absolute_rank_shift = median(detail$absolute_rank_shift),
    topK_dropout_pct           = mean(detail$dropped) * 100
  )
}

## =========================================================
## 3) LOAD + VALIDATE REAL DATA
## =========================================================
if (!file.exists(real_path)) stop("Real data file not found: ", real_path)
real <- read.csv(real_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"class" %in% names(real)) stop("Real data must contain a 'class' column.")

EXPECTED_REAL_ROWS <- 157
EXPECTED_GENE_FEATURES <- 1783
genes <- setdiff(names(real), "class")
if (nrow(real) != EXPECTED_REAL_ROWS) {
  stop("Real: expected ", EXPECTED_REAL_ROWS, " rows, found ", nrow(real), ".")
}
if (length(genes) != EXPECTED_GENE_FEATURES) {
  stop("Real: expected ", EXPECTED_GENE_FEATURES, " genes, found ", length(genes), ".")
}

real$class <- fix_class(real$class, "Real")
real <- force_numeric_finite(real, genes, "Real")

msg("Real data loaded: %d rows x %d genes.", nrow(real), length(genes))
msg("Real class distribution: HC=%d, MDD=%d", sum(real$class == "HC"), sum(real$class == "MDD"))

## =========================================================
## 4) OPTIONAL SYNTHETIC FILES: class balance + observed rankings
## =========================================================
## Defaults below match the archived manuscript files and are used
## only if a synthetic_*.csv file is not found on this machine.
default_balance <- list(
  ARF   = c(hc = 95,  mdd = 105),
  CTGAN = c(hc = 103, mdd = 97),
  TVAE  = c(hc = 103, mdd = 97)
)
synthetic_paths <- list(ARF = arf_path, CTGAN = ctg_path, TVAE = tvae_path)

gen_specs <- list()
for (g in names(synthetic_paths)) {
  p <- synthetic_paths[[g]]
  if (file.exists(p)) {
    syn <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
    syn$class <- fix_class(syn$class, g)
    gen_specs[[g]] <- list(hc = sum(syn$class == "HC"), mdd = sum(syn$class == "MDD"), df = syn)
  } else {
    msg("NOTE: %s synthetic file not found (%s).", g, p)
    msg("      Falling back to the archived manuscript class balance HC=%d, MDD=%d.",
        default_balance[[g]]["hc"], default_balance[[g]]["mdd"])
    gen_specs[[g]] <- list(hc = default_balance[[g]]["hc"], mdd = default_balance[[g]]["mdd"], df = NULL)
  }
}

## =========================================================
## 5) REFERENCE REAL-ONLY RANKING
## =========================================================
msg("\nFitting reference Real-only RF ranking (seed=%d, trees=%d)...", SEED, RANK_TREES)
rank_real_ref <- rf_rank(real, seed_value = SEED, trees = RANK_TREES)

## =========================================================
## 6) OBSERVED Real+Synthetic rankings (computed ONCE, reused below)
## =========================================================
observed_rankings <- list()
for (g in names(gen_specs)) {
  spec <- gen_specs[[g]]
  if (is.null(spec$df)) {
    observed_rankings[[g]] <- NULL
    next
  }
  msg("Fitting OBSERVED Real+%s ranking...", g)
  syn <- spec$df
  syn_genes <- setdiff(names(syn), "class")
  common <- intersect(genes, syn_genes)
  synX <- force_numeric_finite(syn[, c(common, "class")], common, g)
  combined <- bind_rows(real[, c(common, "class")], synX)
  observed_rankings[[g]] <- rf_rank(combined, seed_value = SEED, trees = RANK_TREES)
}

## =========================================================
## SEED BANK (collision-free)
## =========================================================
## Every random draw used by Control A and Control B below pulls its
## seed from ONE global, non-repeating pool (sample.int without
## replacement => every value in the pool is distinct). This replaces
## an earlier hand-rolled additive scheme (SEED + offset*gen_idx +
## 1000*b) that looked independent but was NOT: because the generator
## offset (37000) happened to be an exact multiple of the per-replicate
## step (1000), CTGAN's seed at replicate b was arithmetically IDENTICAL
## to TVAE's seed at replicate (b - 37) for 13 of the 50 replicates,
## silently re-coupling the two "independent" control sets. Drawing
## sequentially from a pre-shuffled, non-repeating pool removes any
## possibility of that kind of arithmetic coincidence, by construction.
## IMPORTANT: draw_seed() must be called in a fixed, deterministic
## order (as written below) for the whole run to be reproducible.
set.seed(SEED + 999000)
n_seeds_needed <- 3 * B_CONTROL_A * length(gen_specs) + B_CONTROL_B + 10
seed_bank <- sample.int(2147483647L, size = n_seeds_needed)
seed_cursor <- 0L
draw_seed <- function() {
  seed_cursor <<- seed_cursor + 1L
  if (seed_cursor > length(seed_bank)) stop("Seed bank exhausted; increase n_seeds_needed.")
  seed_bank[seed_cursor]
}

## =========================================================
## 7) CONTROL A - RESAMPLING CONTROL (no generator)
## =========================================================
msg("\n============================================================")
msg("CONTROL A: resampling 200 rows WITH REPLACEMENT from Real itself")
msg("(class balance matched to each generator), %d replicates for", B_CONTROL_A)
msg("each generator's class-balance condition.")
msg("============================================================")

real_hc  <- real[real$class == "HC", ]
real_mdd <- real[real$class == "MDD", ]

## NOTE: when two generators share the exact same HC/MDD class balance
## (e.g. CTGAN and TVAE both 103/97 in the archived files), their 50
## replicates are still drawn independently (via the seed bank above),
## but both sets are legitimately estimating the SAME underlying
## resampling-noise distribution (the one for a 103/97 class split) -
## they are not expected to be, and are not reported as, a second
## independent form of evidence beyond that one condition.
run_control_a <- function(gen_name, n_hc, n_mdd, B = B_CONTROL_A) {
  stopifnot((n_hc + n_mdd) == N_INJECT)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    if (b %% 10 == 0 || b == 1 || b == B) msg("  %s resampling control: %d/%d", gen_name, b, B)
    set.seed(draw_seed())
    idx_hc  <- sample.int(nrow(real_hc),  size = n_hc,  replace = TRUE)
    set.seed(draw_seed())
    idx_mdd <- sample.int(nrow(real_mdd), size = n_mdd, replace = TRUE)
    resampled <- bind_rows(real_hc[idx_hc, , drop = FALSE], real_mdd[idx_mdd, , drop = FALSE])
    combined  <- bind_rows(real, resampled)
    rank_aug  <- rf_rank(combined, seed_value = draw_seed(), trees = RANK_TREES)
    out[[b]]  <- dilution_vs_real(rank_real_ref, rank_aug) |> mutate(generator = gen_name, replicate = b)
  }
  bind_rows(out)
}

controlA_replicates <- map_dfr(names(gen_specs), function(g) {
  spec <- gen_specs[[g]]
  run_control_a(g, unname(spec$hc), unname(spec$mdd))
})
write_csv(controlA_replicates, file.path(out_dir, "controlA_resampling_replicates.csv"))

## =========================================================
## 8) CONTROL B - RF STOCHASTICITY CONTROL (no data change)
## =========================================================
msg("\n============================================================")
msg("CONTROL B: refit Real-only RF ranking with %d different seeds", B_CONTROL_B)
msg("(same 157 rows every time; isolates pure RF randomness).")
msg("============================================================")

controlB_replicates <- map_dfr(seq_len(B_CONTROL_B), function(b) {
  if (b %% 10 == 0 || b == 1 || b == B_CONTROL_B) msg("  RF-seed control: %d/%d", b, B_CONTROL_B)
  rank_b <- rf_rank(real, seed_value = draw_seed(), trees = RANK_TREES)
  dilution_vs_real(rank_real_ref, rank_b) |> mutate(replicate = b)
})
write_csv(controlB_replicates, file.path(out_dir, "controlB_seed_stability_replicates.csv"))

## =========================================================
## 9) SUPPLEMENTARY CHECK - CORRELATED-GENE SUBSTITUTION
## =========================================================
msg("\n============================================================")
msg("SUPPLEMENTARY CHECK: correlated-gene substitution (Real-only correlations)")
msg("============================================================")

real_gene_matrix <- as.matrix(real[, genes, drop = FALSE])
gene_cor <- cor(real_gene_matrix, method = "spearman")   # Real data only -> non-circular
diag(gene_cor) <- NA

most_correlated_partner <- function(gene) names(which.max(abs(gene_cor[gene, ])))

top_real_genes <- rank_real_ref$gene[seq_len(TOP_K)]

check_substitution <- function(rank_aug, gen_name) {
  dropped_genes <- rank_real_ref |>
    filter(gene %in% top_real_genes) |>
    inner_join(rank_aug |> select(gene, rank_aug = rank), by = "gene") |>
    filter(rank_aug > TOP_K) |>
    pull(gene)
  if (length(dropped_genes) == 0) return(tibble())
  aug_top <- rank_aug$gene[seq_len(TOP_K)]
  partners <- vapply(dropped_genes, most_correlated_partner, character(1))
  tibble(
    generator = gen_name,
    dropped_gene = dropped_genes,
    most_correlated_partner = unname(partners),
    partner_in_augmented_topK = unname(partners) %in% aug_top
  )
}

substitution_results <- map_dfr(names(observed_rankings), function(g) {
  if (is.null(observed_rankings[[g]])) {
    msg("  Skipping %s (synthetic file not found; no observed ranking to check).", g)
    return(NULL)
  }
  check_substitution(observed_rankings[[g]], g)
})

if (nrow(substitution_results) > 0) {
  write_csv(substitution_results, file.path(out_dir, "correlated_gene_substitution.csv"))
  sub_summary <- substitution_results |>
    group_by(generator) |>
    summarise(
      n_dropped_from_topK = n(),
      pct_with_correlated_partner_in_topK = mean(partner_in_augmented_topK) * 100,
      .groups = "drop"
    )
  print(sub_summary)
} else {
  msg("No observed rankings were available; correlated-gene substitution check skipped.")
}

## =========================================================
## 10) SUMMARY TABLE (the one to quote in the response letter)
## =========================================================
## IMPORTANT (two DIFFERENT statistics that must not be read as the
## same quantity):
##   - For the two controls (A and B), "mean_of_replicate_medians" is
##     the MEAN, across B replicates, of each replicate's OWN median
##     absolute rank shift (one median per replicate, then averaged).
##   - For the Observed condition, "median_absolute_rank_shift" is the
##     median computed ONCE, directly, from the single fitted
##     Real+Synthetic ranking (there is only one such ranking per
##     generator, not B replicates of it).
## A control row and the Observed row are therefore never populated in
## the same column; this is intentional, not a missing-data artifact.
summaryA <- controlA_replicates |>
  group_by(generator) |>
  summarise(
    control = "A: resampled Real (no generator)",
    mean_absolute_rank_shift   = mean(mean_absolute_rank_shift),
    mean_of_replicate_medians  = mean(median_absolute_rank_shift),
    median_absolute_rank_shift = NA_real_,
    topK_dropout_pct           = mean(topK_dropout_pct),
    .groups = "drop"
  )

summaryB <- controlB_replicates |>
  summarise(
    generator = "(not applicable)",
    control = "B: RF-seed stochasticity (no data change)",
    mean_absolute_rank_shift   = mean(mean_absolute_rank_shift),
    mean_of_replicate_medians  = mean(median_absolute_rank_shift),
    median_absolute_rank_shift = NA_real_,
    topK_dropout_pct           = mean(topK_dropout_pct)
  )

observed_summary <- map_dfr(names(observed_rankings), function(g) {
  if (is.null(observed_rankings[[g]])) return(NULL)
  dilution_vs_real(rank_real_ref, observed_rankings[[g]]) |>
    mutate(generator = g, control = "Observed Real+Synthetic", mean_of_replicate_medians = NA_real_)
})

final_summary <- bind_rows(summaryA, summaryB, observed_summary) |>
  select(control, generator, mean_absolute_rank_shift, mean_of_replicate_medians,
         median_absolute_rank_shift, topK_dropout_pct)

write_csv(final_summary, file.path(out_dir, "R3_major1_controls_summary.csv"))

msg("\n============================================================")
msg("DONE. Summary table (also saved to R3_major1_controls_summary.csv):")
msg("============================================================")
print(final_summary, n = Inf)

msg("\nAll outputs written under: %s", out_dir)
msg("\nReminder: Control A uses resampled (non-independent) copies of")
msg("existing Real individuals, not new samples; interpret it as the")
msg("effect of RE-WEIGHTING the training set, not as a true new replicate.")
############################################################
# End of R3_rank_displacement_controls.R
############################################################
