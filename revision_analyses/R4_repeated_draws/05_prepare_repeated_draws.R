############################################################
# 05_prepare_repeated_draws.R   (repeated-draws pipeline, step 0)
#
# Creates the seed table for the repeated synthetic-data draws
# (Reviewer 4, point 1; Reviewer 3, Major comment 1).
#
# Design of the whole pipeline
#   * The FULL real data set (gene_with_label.csv, 157 subjects x 1,783
#     genes) is used, exactly as in the main analysis.
#   * Every "draw" repeats the COMPLETE generation procedure from scratch:
#     a new random seed, a new model fit (ARF / CTGAN / TVAE) and a new
#     sample of 200 synthetic subjects. Draws are therefore independent
#     of one another given the 157 real subjects.
#   * Settings are those of the manuscript:
#       ARF   adversarial_rf(100 trees) -> forde(alpha 0.1, epsilon 0.1) -> forge(200)
#       CTGAN 50 epochs, batch 256, PAC 1
#       TVAE  300 epochs, batch 64, embedding 128, (128,128), L2 1e-5, loss factor 2
#     CTGAN_long (optional, expensive) = CTGAN with as many parameter
#     updates as TVAE.
#
# This script only writes <base>/reviewer_revision/R4_repeated_draws/draw_seeds.csv
# (all random seeds of the pipeline, drawn WITHOUT replacement from one
# pool, so no two roles ever share a seed). An existing seed file is
# never overwritten, so the seeds stay fixed across re-runs.
#
# Usage:  Rscript 05_prepare_repeated_draws.R      (or source() it in RStudio)
############################################################
base_dir <- Sys.getenv("RD_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
rd_dir   <- file.path(base_dir, "reviewer_revision", "R4_repeated_draws")
dir.create(rd_dir, recursive = TRUE, showWarnings = FALSE)

GENERATORS <- c("ARF", "CTGAN", "CTGAN_long", "TVAE")
MAX_DRAWS  <- 50L      # seeds are prepared for 50 draws; the default run uses draws 1-30
N_CTRL_B   <- 200L     # seed-only control replicates (real data re-ranked with other RF seeds)
SEED_BASE  <- 5551042

## ---- 1) Check the real data ------------------------------------------------
real_path <- file.path(base_dir, "gene_with_label.csv")
if (!file.exists(real_path)) stop("Real data not found: ", real_path)
real <- read.csv(real_path, check.names = FALSE)
genes <- setdiff(names(real), "class")
if (nrow(real) != 157L || length(genes) != 1783L)
  stop(sprintf("Expected 157 subjects x 1783 genes, found %d x %d.", nrow(real), length(genes)))
if (!setequal(unique(real$class), c(0, 1))) stop("The class column must contain only 0 (HC) and 1 (MDD).")
if (anyNA(real) || any(!is.finite(as.matrix(real[, genes])))) stop("Missing or non-finite values in the real data.")
cat(sprintf("Real data OK: %d subjects (HC %d / MDD %d), %d genes.\n",
            nrow(real), sum(real$class == 0), sum(real$class == 1), length(genes)))

## ---- 2) Seed table ---------------------------------------------------------
f_seeds <- file.path(rd_dir, "draw_seeds.csv")
if (file.exists(f_seeds)) {
  s <- read.csv(f_seeds, stringsAsFactors = FALSE)
  cat(sprintf("draw_seeds.csv already exists (%d seeds) - kept unchanged.\n", nrow(s)))
} else {
  roles <- c("gen", "rank", "ctrlA_hc", "ctrlA_mdd", "ctrlA_rank")
  grid <- do.call(rbind, lapply(GENERATORS, function(g)
    do.call(rbind, lapply(seq_len(MAX_DRAWS), function(d)
      data.frame(generator = g, draw = d, role = roles, stringsAsFactors = FALSE)))))
  ctrlB <- data.frame(generator = "(none)", draw = seq_len(N_CTRL_B), role = "ctrlB_rank", stringsAsFactors = FALSE)
  s <- rbind(grid, ctrlB)
  set.seed(SEED_BASE)
  s$seed <- sample.int(2147483647L, size = nrow(s))      # without replacement => all distinct
  stopifnot(!anyDuplicated(s$seed))
  f_old <- file.path(base_dir, "reviewer_revision", "R3_major2_holdout", "seeds.csv")
  if (file.exists(f_old)) {
    old <- read.csv(f_old)
    if (length(intersect(old$seed, s$seed)) > 0) stop("Seed overlap with the held-out analysis; change SEED_BASE.")
  }
  write.csv(s, f_seeds, row.names = FALSE)
  cat(sprintf("Wrote %d seeds to %s\n", nrow(s), f_seeds))
}
cat("\nNext steps:\n  1) Rscript 05a_draws_arf.R\n  2) python3 05b_draws_ctgan_tvae.py   (long: use caffeinate -i)\n  3) Rscript 05c_evaluate_repeated_draws.R\n")
