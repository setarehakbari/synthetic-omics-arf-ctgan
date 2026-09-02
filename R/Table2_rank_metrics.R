############################################################
## Table2_rank_metrics.R
## Ranking agreement + dilution metrics for Reviewer 1
############################################################

library(tidyverse)
library(ranger)

set.seed(42)

## ---------- Paths ----------
base_dir <- path.expand("~/Desktop/My paper")

real_path <- file.path(base_dir, "gene_with_label.csv")
arf_path  <- file.path(base_dir, "synthetic_rna_arf.csv")

# Your code uses two possible CTGAN filenames, so check both
ctg_candidates <- c(
  file.path(base_dir, "synthetic_ctgan.csv"),
  file.path(base_dir, "synthetic_rna_ctgan.csv")
)

ctg_path <- ctg_candidates[file.exists(ctg_candidates)][1]

if (!file.exists(real_path)) stop("Real data file not found.")
if (!file.exists(arf_path)) stop("ARF data file not found.")
if (is.na(ctg_path)) stop("CTGAN data file not found.")

out_dir <- file.path(base_dir, "reviewer_revision")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Load data ----------
real <- read.csv(real_path, check.names = FALSE)
arf  <- read.csv(arf_path, check.names = FALSE)
ctg  <- read.csv(ctg_path, check.names = FALSE)

## ---------- Class formatting ----------
fix_class <- function(x) {
  z <- as.character(x)
  if (all(unique(z) %in% c("0", "1"))) {
    factor(z, levels = c("0", "1"), labels = c("HC", "MDD"))
  } else {
    factor(z)
  }
}

real$class <- fix_class(real$class)
arf$class  <- fix_class(arf$class)
ctg$class  <- fix_class(ctg$class)

## ---------- Use the SAME genes in all three datasets ----------
genes_real <- setdiff(names(real), "class")
genes_arf  <- setdiff(names(arf), "class")
genes_ctg  <- setdiff(names(ctg), "class")

common <- Reduce(intersect, list(genes_real, genes_arf, genes_ctg))

realX <- real[, c(common, "class")]
arfX  <- arf[, c(common, "class")]
ctgX  <- ctg[, c(common, "class")]

realX[common] <- lapply(realX[common], as.numeric)
arfX[common]  <- lapply(arfX[common], as.numeric)
ctgX[common]  <- lapply(ctgX[common], as.numeric)

## Remove zero-variance genes based on Real data
nzv <- common[sapply(realX[common], sd, na.rm = TRUE) == 0]

if (length(nzv) > 0) {
  common <- setdiff(common, nzv)
  realX <- realX[, c(common, "class")]
  arfX  <- arfX[, c(common, "class")]
  ctgX  <- ctgX[, c(common, "class")]
}

realX <- drop_na(realX)
arfX  <- drop_na(arfX)
ctgX  <- drop_na(ctgX)

cat("Number of genes used:", length(common), "\n")

## ---------- RF ranking ----------
## Use the same RF settings for all three rankings
rf_rank <- function(df, seed_value = 42, trees = 3000) {
  
  m <- ranger(
    class ~ .,
    data = df,
    num.trees = trees,
    importance = "permutation",
    scale.permutation.importance = TRUE,
    seed = seed_value
  )
  
  tibble(
    gene = names(m$variable.importance),
    importance = as.numeric(m$variable.importance)
  ) |>
    arrange(desc(importance), gene) |>
    mutate(rank = row_number())
}

## IMPORTANT:
## Real ranking is calculated ONCE
rank_real <- rf_rank(realX)

rank_realplusARF <- rf_rank(
  bind_rows(realX, arfX)
)

rank_realplusCTGAN <- rf_rank(
  bind_rows(realX, ctgX)
)

## ---------- RBO function ----------
rbo_score <- function(list1, list2, p = 0.95) {
  
  k <- min(length(list1), length(list2))
  
  overlap <- numeric(k)
  
  for (d in seq_len(k)) {
    
    set1 <- list1[1:d]
    set2 <- list2[1:d]
    
    overlap[d] <- length(intersect(set1, set2)) / d
  }
  
  # Extrapolated finite-list RBO
  (1 - p) * sum(overlap * p^(0:(k - 1))) +
    overlap[k] * p^k
}

## ---------- Calculate metrics ----------
calculate_metrics <- function(real_rank,
                              augmented_rank,
                              comparison_name,
                              K = 50,
                              rbo_p = 0.95) {
  
  merged <- real_rank |>
    select(gene, rank_real = rank) |>
    inner_join(
      augmented_rank |>
        select(gene, rank_aug = rank),
      by = "gene"
    )
  
  ## Global ranking agreement
  spearman_rho <- cor(
    merged$rank_real,
    merged$rank_aug,
    method = "spearman"
  )
  
  kendall_tau <- cor(
    merged$rank_real,
    merged$rank_aug,
    method = "kendall"
  )
  
  ## Ranked lists for RBO
  real_list <- real_rank$gene
  aug_list  <- augmented_rank$gene
  
  rbo <- rbo_score(
    real_list,
    aug_list,
    p = rbo_p
  )
  
  ## Top-K sets
  top_real <- head(real_list, K)
  top_aug  <- head(aug_list, K)
  
  jaccard <- length(intersect(top_real, top_aug)) /
    length(union(top_real, top_aug))
  
  ## Dilution only for genes originally in Real Top-K
  detail <- merged |>
    filter(gene %in% top_real) |>
    mutate(
      signed_rank_shift = rank_aug - rank_real,
      absolute_rank_shift = abs(signed_rank_shift),
      dropped = rank_aug > K
    )
  
  result <- tibble(
    Comparison = comparison_name,
    
    Spearman_rho = spearman_rho,
    Kendall_tau = kendall_tau,
    RBO_p_0.95 = rbo,
    Jaccard_K50 = jaccard,
    
    Mean_absolute_rank_shift =
      mean(detail$absolute_rank_shift),
    
    Median_absolute_rank_shift =
      median(detail$absolute_rank_shift),
    
    Top50_dropout_pct =
      mean(detail$dropped) * 100,
    
    ## Secondary directional summaries
    Mean_signed_rank_shift =
      mean(detail$signed_rank_shift),
    
    Median_signed_rank_shift =
      median(detail$signed_rank_shift)
  )
  
  list(summary = result, detail = detail)
}

## ---------- ARF ----------
arf_results <- calculate_metrics(
  rank_real,
  rank_realplusARF,
  "Real vs Real+ARF",
  K = 50,
  rbo_p = 0.95
)

## ---------- CTGAN ----------
ctgan_results <- calculate_metrics(
  rank_real,
  rank_realplusCTGAN,
  "Real vs Real+CTGAN",
  K = 50,
  rbo_p = 0.95
)

## ---------- Final Table ----------
Table2_results <- bind_rows(
  arf_results$summary,
  ctgan_results$summary
)

print(Table2_results)

## Save results
write_csv(
  Table2_results,
  file.path(out_dir, "Table2_ranking_dilution_metrics_K50.csv")
)

write_csv(
  arf_results$detail,
  file.path(out_dir, "dilution_detail_ARF_K50.csv")
)

write_csv(
  ctgan_results$detail,
  file.path(out_dir, "dilution_detail_CTGAN_K50.csv")
)

cat("\nDone.\n")
cat("Results saved in:\n", out_dir, "\n")