############################################################
# 06a_sim_arf.R   (known-signal benchmark, step 1)
# One ARF synthetic data set (200 subjects) per simulated data set, with the
# settings of the manuscript: adversarial_rf(100 trees), forde(alpha 0.1,
# epsilon 0.1), forge(200). Finished sims are skipped (restartable).
# Usage: Rscript 06a_sim_arf.R        (KS_QUICK=1: smoke test in a separate folder)
# Output: draws/ARF/sim_XX.csv (+ log_XX.csv). Simulated data: no privacy restrictions.
############################################################
ks_dir <- Sys.getenv("KS_DIR", unset = ".")
QUICK <- identical(Sys.getenv("KS_QUICK", "0"), "1")
FIRST <- as.integer(Sys.getenv("KS_FIRST", "1")); LAST <- as.integer(Sys.getenv("KS_LAST", "10"))
N_SYNTH <- 200L; GEN <- "ARF"
gen_dir <- file.path(ks_dir, if (QUICK) "draws_QUICK" else "draws", GEN); dir.create(gen_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(arf))
seeds <- read.csv(file.path(ks_dir, "seeds.csv"), stringsAsFactors = FALSE)
for (i in FIRST:LAST) {
  out <- file.path(gen_dir, sprintf("sim_%02d.csv", i))
  if (file.exists(out)) { cat(sprintf("ARF sim %02d: exists, skipping\n", i)); next }
  real <- read.csv(file.path(ks_dir, "sims", sprintf("sim_%02d.csv", i)), check.names = FALSE)
  genes <- setdiff(names(real), "class"); stopifnot(nrow(real) == 157, setequal(unique(real$class), c(0, 1)))
  real$class <- factor(as.character(real$class), levels = c("0", "1"))
  set.seed(seeds$seed[seeds$sim == i & seeds$role == "gen_ARF"])
  t0 <- proc.time()[["elapsed"]]
  m   <- adversarial_rf(real, num_trees = if (QUICK) 20L else 100L)
  psi <- forde(m, real, alpha = 0.1, epsilon = 0.1)
  syn <- as.data.frame(forge(psi, n_synth = N_SYNTH), check.names = FALSE)
  secs <- proc.time()[["elapsed"]] - t0
  stopifnot(nrow(syn) == N_SYNTH, setequal(names(syn), names(real))); syn <- syn[, names(real), drop = FALSE]
  syn$class <- suppressWarnings(as.integer(as.character(syn$class)))
  for (g in genes) syn[[g]] <- suppressWarnings(as.numeric(syn[[g]]))
  if (anyNA(syn) || !all(syn$class %in% c(0L, 1L)) || any(!is.finite(as.matrix(syn[, genes])))) stop("Invalid synthetic data, sim ", i)
  tmp <- paste0(out, ".part"); write.csv(syn, tmp, row.names = FALSE); file.rename(tmp, out)
  write.csv(data.frame(sim = i, generator = GEN, n_synth = nrow(syn), synth_class0 = sum(syn$class == 0), synth_class1 = sum(syn$class == 1),
                       fit_seconds = round(secs, 1), arf_version = as.character(packageVersion("arf")), quick = QUICK),
            file.path(gen_dir, sprintf("log_%02d.csv", i)), row.names = FALSE)
  cat(sprintf("ARF sim %02d done in %.0f s (class 0/1 = %d/%d)\n", i, secs, sum(syn$class == 0), sum(syn$class == 1)))
}
cat("ARF finished.\n")
