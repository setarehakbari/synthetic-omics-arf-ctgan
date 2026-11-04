############################################################
# 06d_sim_sensitivity_shift_only.R   (known-signal benchmark, sensitivity analysis)
#
# Top-50 dropout is an exact function of Top-50 Jaccard similarity
# (Jaccard = (50 - d) / (50 + d) for d dropped genes), so the "dilution" and the
# "agreement" sets of 06c_sim_evaluate.R share one measure. This script repeats the
# incremental R^2 comparison with only the mean and median absolute rank shift as
# dilution summaries, and reports the R^2 of each single summary.
#
# Usage:  Rscript 06d_sim_sensitivity_shift_only.R      (reads results/per_ranking_metrics.csv)
# Output: results/table_incremental_R2_shift_only.csv, results/table_single_summary_R2.csv
############################################################
ks_dir <- Sys.getenv("KS_DIR", unset = ".")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
res <- read.csv(file.path(ks_dir, "results", "per_ranking_metrics.csv"), check.names = FALSE)
nonref <- res %>% filter(condition != "Reference")
K <- 50; dd <- nonref$dropout_pct / 100 * K
stopifnot(max(abs(nonref$jaccard_top50 - (K - dd) / (K + dd)), na.rm = TRUE) < 1e-8)   # exact relation
set.seed(20261029); NB <- 2000
r2 <- function(df, outcome, preds) {
  d <- df %>% group_by(sim) %>% mutate(across(all_of(c(outcome, preds)), ~ .x - mean(.x))) %>% ungroup()   # simulation fixed effects
  X <- as.matrix(d[, preds]); y <- d[[outcome]]; 1 - sum(residuals(lm(y ~ X - 1))^2) / sum(y^2)
}
SETS <- list(rankagree = c("spearman", "rbo_p0.95", "jaccard_top50"), shift = c("mean_abs_shift", "median_abs_shift"),
             both = c("spearman", "rbo_p0.95", "jaccard_top50", "mean_abs_shift", "median_abs_shift"))
single <- c("mean_abs_shift", "median_abs_shift", "dropout_pct", "spearman", "rbo_p0.95", "jaccard_top50")
run <- function(df, outcome, label) {
  one <- function(d) c(sapply(SETS, function(p) r2(d, outcome, p)), sapply(single, function(m) r2(d, outcome, m)))
  est <- one(df); su <- unique(df$sim)
  bs <- t(replicate(NB, { s <- sample(su, replace = TRUE); one(bind_rows(lapply(seq_along(s), function(j) mutate(df[df$sim == s[j], ], sim = j)))) }))
  ci <- function(v) quantile(v, c(.025, .975))
  list(inc = tibble(subset = label, outcome = outcome, R2_agree = est["rankagree"], R2_shift = est["shift"], R2_both = est["both"],
                    add_shift_to_agree = est["both"] - est["rankagree"], lo1 = ci(bs[, "both"] - bs[, "rankagree"])[1], hi1 = ci(bs[, "both"] - bs[, "rankagree"])[2],
                    add_agree_to_shift = est["both"] - est["shift"], lo2 = ci(bs[, "both"] - bs[, "shift"])[1], hi2 = ci(bs[, "both"] - bs[, "shift"])[2]),
       single = tibble(subset = label, outcome = outcome, measure = single, R2 = est[single],
                       lo = sapply(single, function(m) ci(bs[, m])[1]), hi = sapply(single, function(m) ci(bs[, m])[2])))
}
nc <- nonref %>% filter(condition != "Control C")
outs <- list(run(nonref, "loss_tp", "all rankings"), run(nc, "loss_tp", "without Control C"),
             run(nonref, "loss_auroc", "all rankings"), run(nc, "loss_auroc", "without Control C"))
write.csv(bind_rows(lapply(outs, `[[`, "inc")), file.path(ks_dir, "results", "table_incremental_R2_shift_only.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(outs, `[[`, "single")), file.path(ks_dir, "results", "table_single_summary_R2.csv"), row.names = FALSE)
print(as.data.frame(bind_rows(lapply(outs, `[[`, "inc"))), digits = 3)
cat("Results written to:", file.path(ks_dir, "results"), "\n")
