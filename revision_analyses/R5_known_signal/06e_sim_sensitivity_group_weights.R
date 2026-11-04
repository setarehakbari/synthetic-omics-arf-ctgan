############################################################
# 06e_sim_sensitivity_group_weights.R   (known-signal benchmark, sensitivity analysis)
#
# In 06c_sim_evaluate.R every ranking is one row, so Control B (30 re-rankings per
# simulation) contributes about 77% of the non-reference rows. This script repeats
# the incremental R^2 comparison giving every condition x generator group within a
# simulation the same total weight (each ranking is weighted 1 / group size;
# Observed: one ranking per generator, Control A: 3 per generator, Control B: 30, Control C: 3 per simulation),
# using weighted within-simulation least squares.
#
# Usage:  Rscript 06e_sim_sensitivity_group_weights.R      (reads results/per_ranking_metrics.csv)
# Output: results/table_incremental_R2_group_weights.csv
############################################################
ks_dir <- Sys.getenv("KS_DIR", unset = ".")
suppressPackageStartupMessages({ library(dplyr) })
res <- read.csv(file.path(ks_dir, "results", "per_ranking_metrics.csv"), check.names = FALSE)
nonref <- res %>% filter(condition != "Reference") %>%
  group_by(sim, condition, generator) %>% mutate(w = 1 / n()) %>% ungroup()
set.seed(20261029); NB <- 2000
INC_AGR <- c("spearman", "rbo_p0.95", "jaccard_top50")
INC_DISP <- c("mean_abs_shift", "median_abs_shift", "dropout_pct"); SHIFT <- c("mean_abs_shift", "median_abs_shift")
wr2 <- function(df, outcome, preds) {
  d <- df %>% group_by(sim) %>% mutate(across(all_of(c(outcome, preds)), ~ .x - weighted.mean(.x, w))) %>% ungroup()
  X <- as.matrix(d[, preds]); y <- d[[outcome]]; w <- d$w
  fit <- lm.wfit(X, y, w); 1 - sum(w * fit$residuals^2) / sum(w * y^2)
}
run <- function(df, outcome, label) {
  one <- function(d) c(agree = wr2(d, outcome, INC_AGR), disp = wr2(d, outcome, INC_DISP), shift = wr2(d, outcome, SHIFT),
                       both = wr2(d, outcome, c(INC_AGR, INC_DISP)), both_shift = wr2(d, outcome, c(INC_AGR, SHIFT)))
  est <- one(df); su <- unique(df$sim)
  bs <- t(replicate(NB, { s <- sample(su, replace = TRUE); one(bind_rows(lapply(seq_along(s), function(j) mutate(df[df$sim == s[j], ], sim = j)))) }))
  ci <- function(v) quantile(v, c(.025, .975))
  tibble(subset = label, outcome = outcome, R2_agree = est[["agree"]], R2_dilution = est[["disp"]], R2_both = est[["both"]],
         add_dilution_to_agree = est[["both"]] - est[["agree"]], lo1 = ci(bs[, "both"] - bs[, "agree"])[[1]], hi1 = ci(bs[, "both"] - bs[, "agree"])[[2]],
         add_agree_to_dilution = est[["both"]] - est[["disp"]], lo2 = ci(bs[, "both"] - bs[, "disp"])[[1]], hi2 = ci(bs[, "both"] - bs[, "disp"])[[2]],
         add_shift_to_agree = est[["both_shift"]] - est[["agree"]], lo3 = ci(bs[, "both_shift"] - bs[, "agree"])[[1]], hi3 = ci(bs[, "both_shift"] - bs[, "agree"])[[2]],
         add_agree_to_shift = est[["both_shift"]] - est[["shift"]], lo4 = ci(bs[, "both_shift"] - bs[, "shift"])[[1]], hi4 = ci(bs[, "both_shift"] - bs[, "shift"])[[2]])
}
nc <- nonref %>% filter(condition != "Control C")
out <- bind_rows(run(nonref, "loss_tp", "all rankings (incl. Control C)"), run(nc, "loss_tp", "without Control C"))
write.csv(out, file.path(ks_dir, "results", "table_incremental_R2_group_weights.csv"), row.names = FALSE)
print(as.data.frame(out), digits = 3)
cat("Results written to:", file.path(ks_dir, "results"), "\n")
