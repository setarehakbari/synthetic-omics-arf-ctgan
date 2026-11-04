############################################################
# 05a_draws_arf.R   (repeated-draws pipeline, step 1)
#
# Generates repeated, independent ARF draws from the FULL real data.
# One draw = new seed + new adversarial_rf fit + forde + forge(200).
# Settings as in the manuscript (adversarial_rf 100 trees,
# forde alpha = 0.1, epsilon = 0.1, 200 synthetic subjects).
#
# Usage
#   Rscript 05a_draws_arf.R                   # draws 1-30
#   RD_FIRST=1 RD_LAST=5 Rscript 05a_draws_arf.R
#   RD_QUICK=1 RD_LAST=2 Rscript 05a_draws_arf.R   # smoke test, separate folder
# Finished draws are skipped, so the script can be stopped and restarted.
#
# Output: <base>/reviewer_revision/R4_repeated_draws/draws/ARF/draw_XX.csv (+ log_XX.csv)
# These files contain synthetic profiles derived from patient data:
# keep them OUT of any public repository.
############################################################
base_dir  <- Sys.getenv("RD_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
rd_dir    <- file.path(base_dir, "reviewer_revision", "R4_repeated_draws")
QUICK     <- identical(Sys.getenv("RD_QUICK", "0"), "1")
OVERWRITE <- identical(Sys.getenv("RD_OVERWRITE", "0"), "1")
FIRST <- as.integer(Sys.getenv("RD_FIRST", "1")); LAST <- as.integer(Sys.getenv("RD_LAST", "30"))
N_SYNTH <- 200L; ARF_NUM_TREES <- 100L; FORDE_ALPHA <- 0.1; FORDE_EPSILON <- 0.1
GEN <- "ARF"
gen_dir <- file.path(rd_dir, if (QUICK) "draws_QUICK" else "draws", GEN)
dir.create(gen_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("arf", quietly = TRUE)) stop("Package 'arf' is required: install.packages('arf')")
suppressPackageStartupMessages(library(arf))

seeds <- read.csv(file.path(rd_dir, "draw_seeds.csv"), stringsAsFactors = FALSE)
get_seed <- function(d) {
  s <- seeds$seed[seeds$generator == GEN & seeds$draw == d & seeds$role == "gen"]
  if (length(s) != 1) stop("No seed for ", GEN, " draw ", d, ". Run 05_prepare_repeated_draws.R first.")
  s
}

real <- read.csv(file.path(base_dir, "gene_with_label.csv"), check.names = FALSE)
genes <- setdiff(names(real), "class")
stopifnot(nrow(real) == 157, length(genes) == 1783, all(is.finite(as.matrix(real[, genes]))),
          setequal(unique(real$class), c(0, 1)))
real$class <- factor(as.character(real$class), levels = c("0", "1"))

for (d in FIRST:LAST) {
  out <- file.path(gen_dir, sprintf("draw_%02d.csv", d))
  if (file.exists(out) && !OVERWRITE) { cat(sprintf("ARF draw %02d: exists, skipping\n", d)); next }
  set.seed(get_seed(d))
  t0 <- proc.time()[["elapsed"]]
  arf_model <- adversarial_rf(real, num_trees = ARF_NUM_TREES)
  psi       <- forde(arf_model, real, alpha = FORDE_ALPHA, epsilon = FORDE_EPSILON)
  syn       <- as.data.frame(forge(psi, n_synth = N_SYNTH), check.names = FALSE)
  secs <- proc.time()[["elapsed"]] - t0

  stopifnot(nrow(syn) == N_SYNTH, setequal(names(syn), names(real)))
  syn <- syn[, names(real), drop = FALSE]
  syn$class <- suppressWarnings(as.integer(as.character(syn$class)))
  for (g in genes) syn[[g]] <- suppressWarnings(as.numeric(syn[[g]]))
  if (anyNA(syn) || !all(syn$class %in% c(0L, 1L)) || any(!is.finite(as.matrix(syn[, genes]))))
    stop("Invalid synthetic data in ARF draw ", d)
  tmp <- paste0(out, ".part"); write.csv(syn, tmp, row.names = FALSE); file.rename(tmp, out)   # no half-written files
  write.csv(data.frame(
    draw = d, generator = GEN, seed_role = "gen", n_real = nrow(real), n_synth = nrow(syn),
    synth_class0 = sum(syn$class == 0), synth_class1 = sum(syn$class == 1),
    arf_trees = ARF_NUM_TREES, forde_alpha = FORDE_ALPHA, forde_epsilon = FORDE_EPSILON,
    fit_seconds = round(secs, 1), arf_version = as.character(packageVersion("arf")),
    quick = QUICK, stringsAsFactors = FALSE),
    file.path(gen_dir, sprintf("log_%02d.csv", d)), row.names = FALSE)
  cat(sprintf("ARF draw %02d done in %.0f s (synthetic class 0/1 = %d/%d)\n",
              d, secs, sum(syn$class == 0), sum(syn$class == 1)))
}
cat("ARF finished.\n")
