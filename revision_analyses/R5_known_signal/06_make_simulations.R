############################################################
# 06_make_simulations.R   (known-signal benchmark, step 0)
#
# Simulates S data sets with KNOWN informative genes using the npdro
# simulator (createSimulation2, settings suggested by B. McKinney).
# No patient data are used. The simulated data sets mimic the dimensions
# of the real study (157 subjects, ~1,800 genes, 79 / 78 per class).
#
#   * 5% of the genes are signals: 45 "main effect" genes (class mean shift)
#     and 45 "interaction" genes (effect only through gene-gene interactions;
#     they are rarely recovered by random-forest permutation importance).
#   * main.bias = interaction.bias = 1.5, hi.cor = 0.95, lo.cor = 0.2,
#     mix.type = "main-interactionScalefree", pct.mixed = 0.5.
#
# The ready-made output is part of the package that you received
# (sims/, truth.csv, sim_info.csv, seeds.csv); re-run this script only to
# regenerate it (needs: remotes::install_github("insilico/npdro")).
#
# Usage:  Rscript 06_make_simulations.R
# Output (never overwritten): sims/sim_XX.csv, truth.csv, sim_info.csv, seeds.csv
############################################################
ks_dir <- Sys.getenv("KS_DIR", unset = ".")
dir.create(file.path(ks_dir, "sims"), recursive = TRUE, showWarnings = FALSE)
N_SIM  <- as.integer(Sys.getenv("KS_N_SIM", "10"))
N_VAR  <- 1800L; N_SAMP <- 157L; N_HC <- 79L
N_CTRL_B <- 30L; N_CTRL_C <- 3L           # controls per simulated data set
GENS   <- c("ARF", "CTGAN", "TVAE")
SEED_BASE <- 7770042
suppressPackageStartupMessages(library(npdro))

## ---- seed table (all seeds distinct, drawn without replacement) --------------------------
f_seeds <- file.path(ks_dir, "seeds.csv")
if (!file.exists(f_seeds)) {
  roles <- c("sim", "rank_ref",
             paste0("gen_", GENS), paste0("rank_obs_", GENS),
             as.vector(sapply(GENS, function(g) paste0(c("ctrlA_hc_", "ctrlA_mdd_", "ctrlA_rank_"), g))),
             paste0("ctrlB_rank_", seq_len(N_CTRL_B)),
             as.vector(sapply(seq_len(N_CTRL_C), function(k) paste0(c("ctrlC_idx_", "ctrlC_perm_", "ctrlC_rank_"), k))))
  s <- expand.grid(sim = seq_len(N_SIM), role = roles, stringsAsFactors = FALSE)
  set.seed(SEED_BASE); s$seed <- sample.int(2147483647L, nrow(s)); stopifnot(!anyDuplicated(s$seed))
  write.csv(s, f_seeds, row.names = FALSE); cat("Wrote", nrow(s), "seeds\n")
}
seeds <- read.csv(f_seeds, stringsAsFactors = FALSE)
get_seed <- function(i, role) seeds$seed[seeds$sim == i & seeds$role == role]

## ---- simulate -----------------------------------------------------------------------------
truth <- list(); info <- list()
for (i in seq_len(N_SIM)) {
  out <- file.path(ks_dir, "sims", sprintf("sim_%02d.csv", i))
  if (file.exists(out)) { cat(sprintf("sim %02d exists, skipping\n", i)); next }
  set.seed(get_seed(i, "sim"))
  ds <- createSimulation2(num.samples = N_SAMP, num.variables = N_VAR, pct.imbalance = N_HC / N_SAMP,
                          pct.signals = 0.05, main.bias = 1.5, interaction.bias = 1.5, hi.cor = 0.95, lo.cor = 0.2,
                          mix.type = "main-interactionScalefree", label = "class", sim.type = "mixed",
                          graph.structure = NULL, pct.mixed = 0.5, pct.train = 1, pct.holdout = 0,
                          pct.validation = 0, plot.graph = FALSE, verbose = FALSE)
  d <- ds$train
  d$class <- ifelse(d$class == 1, 1L, 0L)
  genes <- setdiff(names(d), "class")
  stopifnot(nrow(d) == N_SAMP, all(is.finite(as.matrix(d[, genes]))), length(ds$signal.names) == 90)
  tmp <- paste0(out, ".part"); write.csv(d, tmp, row.names = FALSE); file.rename(tmp, out)
  tr <- data.frame(sim = i, gene = ds$signal.names,
                   type = ifelse(grepl("^main", ds$signal.names), "main", "interaction"), stringsAsFactors = FALSE)
  write.csv(tr, file.path(ks_dir, "sims", sprintf("truth_%02d.csv", i)), row.names = FALSE)
  cat(sprintf("sim %02d: %d x %d, class 0/1 = %d/%d, signals = %d (%d main, %d interaction)\n", i, nrow(d), length(genes),
              sum(d$class == 0), sum(d$class == 1), nrow(tr), sum(tr$type == "main"), sum(tr$type == "interaction")))
}
tf <- list.files(file.path(ks_dir, "sims"), "^truth_\\d+\\.csv$", full.names = TRUE)
write.csv(do.call(rbind, lapply(tf, read.csv)), file.path(ks_dir, "truth.csv"), row.names = FALSE)
cat("Done. Next: Rscript 06a_sim_arf.R ; python3 06b_sim_ctgan_tvae.py ; Rscript 06c_sim_evaluate.R\n")
