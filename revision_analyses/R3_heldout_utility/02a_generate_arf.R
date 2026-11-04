############################################################
# 02a_generate_arf.R   (stage 2a of 3)
#
# Fits an Adversarial Random Forest (ARF) on the TRAINING subjects
# of each split only, and generates 200 synthetic subjects.
# Settings follow the manuscript (ARF_generate.R):
#   adversarial_rf(num_trees = 100) -> forde(alpha = 0.1,
#   epsilon = 0.1) -> forge(n_synth = 200).
#
# Usage
#   Rscript 02a_generate_arf.R               # all splits
#   HOLDOUT_SPLITS=1,2 Rscript 02a_generate_arf.R   # selected splits
# Existing outputs are kept unless HOLDOUT_OVERWRITE=1.
#
# Output per split: splits/split_XX/synthetic_ARF.csv and
#                   splits/split_XX/log_ARF.csv
############################################################
base_dir <- Sys.getenv("HOLDOUT_BASE_DIR", unset = Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = "."))
out_dir  <- file.path(base_dir, "reviewer_revision", "R3_major2_holdout")
N_SYNTH <- 200L; ARF_NUM_TREES <- 100L; FORDE_ALPHA <- 0.1; FORDE_EPSILON <- 0.1
OVERWRITE <- identical(Sys.getenv("HOLDOUT_OVERWRITE", "0"), "1")

if (!requireNamespace("arf", quietly = TRUE)) stop("Package 'arf' is required.")
suppressPackageStartupMessages(library(arf))

seeds <- read.csv(file.path(out_dir, "seeds.csv"), stringsAsFactors = FALSE)
get_seed <- function(i, role) seeds$seed[seeds$split == i & seeds$role == role]

all_splits <- sort(as.integer(sub("split_", "", list.files(file.path(out_dir, "splits"), pattern = "^split_"))))
sel <- Sys.getenv("HOLDOUT_SPLITS", unset = "")
splits <- if (nzchar(sel)) as.integer(strsplit(sel, ",")[[1]]) else all_splits

for (i in splits) {
  d   <- file.path(out_dir, "splits", sprintf("split_%02d", i))
  out <- file.path(d, "synthetic_ARF.csv")
  if (file.exists(out) && !OVERWRITE) { cat(sprintf("split %02d: ARF output exists, skipping\n", i)); next }

  train <- read.csv(file.path(d, "train.csv"), check.names = FALSE)
  genes <- setdiff(names(train), "class")
  stopifnot(all(is.finite(as.matrix(train[, genes]))), setequal(unique(train$class), c(0, 1)))
  train$class <- factor(as.character(train$class), levels = c("0", "1"))

  set.seed(get_seed(i, "gen_ARF"))
  t0 <- proc.time()[["elapsed"]]
  arf_model <- adversarial_rf(train, num_trees = ARF_NUM_TREES)
  psi       <- forde(arf_model, train, alpha = FORDE_ALPHA, epsilon = FORDE_EPSILON)
  syn       <- as.data.frame(forge(psi, n_synth = N_SYNTH), check.names = FALSE)
  secs <- proc.time()[["elapsed"]] - t0

  stopifnot(nrow(syn) == N_SYNTH, setequal(names(syn), names(train)))
  syn <- syn[, names(train), drop = FALSE]
  syn$class <- suppressWarnings(as.integer(as.character(syn$class)))
  for (g in genes) syn[[g]] <- suppressWarnings(as.numeric(syn[[g]]))
  if (anyNA(syn) || !all(syn$class %in% c(0L, 1L)) || any(!is.finite(as.matrix(syn[, genes]))))
    stop("Invalid synthetic data in split ", i)
  write.csv(syn, out, row.names = FALSE)

  write.csv(data.frame(
    split = i, generator = "ARF", n_train = nrow(train), n_synth = nrow(syn),
    synth_class0 = sum(syn$class == 0), synth_class1 = sum(syn$class == 1),
    arf_trees = ARF_NUM_TREES, forde_alpha = FORDE_ALPHA, forde_epsilon = FORDE_EPSILON,
    fit_seconds = round(secs, 1), arf_version = as.character(packageVersion("arf")),
    stringsAsFactors = FALSE), file.path(d, "log_ARF.csv"), row.names = FALSE)
  cat(sprintf("split %02d: ARF done in %.0f s (synthetic class 0/1 = %d/%d)\n",
              i, secs, sum(syn$class == 0), sum(syn$class == 1)))
}
