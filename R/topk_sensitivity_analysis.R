############################################################
# topk_sensitivity_analysis.R
# Sensitivity analysis for Reviewer Minor Comment 5
# Uses SAVED full RF and NPDR ranking files only.
# No RF/NPDR model refitting and no bootstrap rerun.
############################################################

suppressPackageStartupMessages(library(tidyverse))

K_VALUES <- c(25, 50, 100)
GENERATORS <- c("ARF", "CTGAN", "TVAE")

out_dir <- path.expand("~/Desktop/My paper/reviewer_revision/all_generators")
details_dir <- file.path(out_dir, "details")
tables_dir  <- file.path(out_dir, "tables")

if (!dir.exists(details_dir)) stop("details/ folder not found: ", details_dir)
if (!dir.exists(tables_dir))  stop("tables/ folder not found: ", tables_dir)

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

read_rank <- function(path) {
  if (!file.exists(path)) stop("Missing ranking file: ", path)
  d <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("gene", "rank")
  miss <- setdiff(required, names(d))
  if (length(miss) > 0) stop("File ", basename(path), " is missing: ", paste(miss, collapse = ", "))
  d <- d |>
    transmute(gene = as.character(gene), rank = as.numeric(rank)) |>
    arrange(rank, gene)
  if (anyDuplicated(d$gene)) stop("Duplicated genes in: ", basename(path))
  if (anyNA(d$rank)) stop("Missing ranks in: ", basename(path))
  d
}

check_same_gene_set <- function(reference, x, label) {
  if (!setequal(reference$gene, x$gene)) stop("Gene set mismatch in ", label)
}

alignment_rows <- function(method, generator, real_rank, syn_rank, comb_rank, K) {
  top_real <- head(real_rank$gene, K)
  top_syn  <- head(syn_rank$gene, K)
  top_comb <- head(comb_rank$gene, K)
  tibble(
    analysis = paste0(method, " feature alignment"),
    method = method,
    generator = generator,
    comparison = c("Real vs Synthetic", "Real vs Real+Synthetic", "Synthetic vs Real+Synthetic"),
    K = K,
    n_overlap = c(length(intersect(top_real, top_syn)),
                  length(intersect(top_real, top_comb)),
                  length(intersect(top_syn, top_comb))),
    jaccard = c(jaccard(top_real, top_syn),
                jaccard(top_real, top_comb),
                jaccard(top_syn, top_comb)),
    mean_absolute_rank_shift = NA_real_,
    median_absolute_rank_shift = NA_real_,
    topK_dropout_pct = NA_real_,
    mean_signed_rank_shift = NA_real_,
    median_signed_rank_shift = NA_real_
  )
}

dilution_row <- function(generator, real_rank, aug_rank, K) {
  merged <- real_rank |>
    select(gene, rank_real = rank) |>
    inner_join(aug_rank |> select(gene, rank_aug = rank), by = "gene")
  if (nrow(merged) != nrow(real_rank)) stop("Ranking merge failed for Real vs Real+", generator)

  top_real <- head(real_rank$gene, K)
  top_aug  <- head(aug_rank$gene, K)
  detail <- merged |>
    filter(gene %in% top_real) |>
    mutate(signed_rank_shift = rank_aug - rank_real,
           absolute_rank_shift = abs(signed_rank_shift),
           dropped = rank_aug > K)

  tibble(
    analysis = "Dilution",
    method = "RF permutation importance",
    generator = generator,
    comparison = paste0("Real vs Real+", generator),
    K = K,
    n_overlap = length(intersect(top_real, top_aug)),
    jaccard = jaccard(top_real, top_aug),
    mean_absolute_rank_shift = mean(detail$absolute_rank_shift),
    median_absolute_rank_shift = median(detail$absolute_rank_shift),
    topK_dropout_pct = mean(detail$dropped) * 100,
    mean_signed_rank_shift = mean(detail$signed_rank_shift),
    median_signed_rank_shift = median(detail$signed_rank_shift)
  )
}

assert_close <- function(x, y, label, tol = 1e-10) {
  if (length(x) != length(y) || any(abs(x - y) > tol, na.rm = TRUE)) {
    stop("K=50 validation FAILED for ", label,
         ". Sensitivity calculations do not reproduce the primary analysis.")
  }
}

## Load saved RF rankings
rf_real <- read_rank(file.path(details_dir, "05_rf_importance_Real.csv"))
rf_syn <- setNames(vector("list", length(GENERATORS)), GENERATORS)
rf_comb <- setNames(vector("list", length(GENERATORS)), GENERATORS)
for (g in GENERATORS) {
  rf_syn[[g]]  <- read_rank(file.path(details_dir, paste0("05_rf_importance_", g, ".csv")))
  rf_comb[[g]] <- read_rank(file.path(details_dir, paste0("05_rf_importance_RealPlus", g, ".csv")))
  check_same_gene_set(rf_real, rf_syn[[g]], paste0("RF ", g))
  check_same_gene_set(rf_real, rf_comb[[g]], paste0("RF Real+", g))
}

## Load saved NPDR rankings
npdr_real <- read_rank(file.path(details_dir, "06_npdr_importance_Real.csv"))
npdr_syn <- setNames(vector("list", length(GENERATORS)), GENERATORS)
npdr_comb <- setNames(vector("list", length(GENERATORS)), GENERATORS)
for (g in GENERATORS) {
  npdr_syn[[g]]  <- read_rank(file.path(details_dir, paste0("06_npdr_importance_", g, ".csv")))
  npdr_comb[[g]] <- read_rank(file.path(details_dir, paste0("06_npdr_importance_RealPlus", g, ".csv")))
  check_same_gene_set(npdr_real, npdr_syn[[g]], paste0("NPDR ", g))
  check_same_gene_set(npdr_real, npdr_comb[[g]], paste0("NPDR Real+", g))
}

## Compute K = 25, 50, 100
rf_results <- bind_rows(lapply(K_VALUES, function(K) {
  bind_rows(lapply(GENERATORS, function(g) {
    alignment_rows("RF", g, rf_real, rf_syn[[g]], rf_comb[[g]], K)
  }))
}))

npdr_results <- bind_rows(lapply(K_VALUES, function(K) {
  bind_rows(lapply(GENERATORS, function(g) {
    alignment_rows("NPDR", g, npdr_real, npdr_syn[[g]], npdr_comb[[g]], K)
  }))
}))

dilution_results <- bind_rows(lapply(K_VALUES, function(K) {
  bind_rows(lapply(GENERATORS, function(g) {
    dilution_row(g, rf_real, rf_comb[[g]], K)
  }))
}))

## Validate K=50 against the primary manuscript tables
primary_rf_path   <- file.path(tables_dir, "05_rf_feature_alignment_summary.csv")
primary_npdr_path <- file.path(tables_dir, "06_npdr_feature_alignment_summary.csv")
primary_dil_path  <- file.path(tables_dir, "07_ranking_dilution_summary_K50.csv")
if (!all(file.exists(c(primary_rf_path, primary_npdr_path, primary_dil_path)))) {
  stop("One or more primary K=50 summary tables are missing; cannot validate sensitivity scope.")
}

primary_rf <- readr::read_csv(primary_rf_path, show_col_types = FALSE) |> arrange(generator, comparison)
new_rf50 <- rf_results |> filter(K == 50) |> arrange(generator, comparison)
assert_close(new_rf50$jaccard, primary_rf$jaccard, "RF Jaccard")
assert_close(new_rf50$n_overlap, primary_rf$n_overlap, "RF overlap")

primary_npdr <- readr::read_csv(primary_npdr_path, show_col_types = FALSE) |> arrange(generator, comparison)
new_npdr50 <- npdr_results |> filter(K == 50) |> arrange(generator, comparison)
assert_close(new_npdr50$jaccard, primary_npdr$jaccard, "NPDR Jaccard")
assert_close(new_npdr50$n_overlap, primary_npdr$n_overlap, "NPDR overlap")

primary_dil <- readr::read_csv(primary_dil_path, show_col_types = FALSE) |>
  mutate(generator = sub("^Real vs Real\\+", "", Comparison)) |>
  arrange(generator)
new_dil50 <- dilution_results |> filter(K == 50) |> arrange(generator)
assert_close(new_dil50$jaccard, primary_dil$Jaccard_K50, "Dilution Jaccard")
assert_close(new_dil50$mean_absolute_rank_shift, primary_dil$Mean_absolute_rank_shift,
             "Dilution mean absolute rank shift")
assert_close(new_dil50$median_absolute_rank_shift, primary_dil$Median_absolute_rank_shift,
             "Dilution median absolute rank shift")
assert_close(new_dil50$topK_dropout_pct, primary_dil$Top50_dropout_pct,
             "Dilution Top-50 dropout")

message("K=50 VALIDATION PASSED: sensitivity calculations exactly reproduce the primary analysis.")

## Save outputs
readr::write_csv(rf_results, file.path(tables_dir, "10a_RF_TopK_sensitivity_K25_K50_K100.csv"))
readr::write_csv(npdr_results, file.path(tables_dir, "10b_NPDR_TopK_sensitivity_K25_K50_K100.csv"))
readr::write_csv(dilution_results, file.path(tables_dir, "10c_Dilution_TopK_sensitivity_K25_K50_K100.csv"))

all_results <- bind_rows(rf_results, npdr_results, dilution_results)
readr::write_csv(all_results, file.path(tables_dir, "Supplementary_Table_S1_TopK_sensitivity.csv"))

compact <- bind_rows(
  rf_results |>
    filter(comparison == "Real vs Real+Synthetic") |>
    select(analysis, generator, K, n_overlap, jaccard),
  npdr_results |>
    filter(comparison == "Real vs Real+Synthetic") |>
    select(analysis, generator, K, n_overlap, jaccard),
  dilution_results |>
    select(analysis, generator, K, n_overlap, jaccard,
           mean_absolute_rank_shift, median_absolute_rank_shift, topK_dropout_pct)
)
readr::write_csv(compact, file.path(tables_dir, "Supplementary_Table_S1_TopK_sensitivity_compact.csv"))

cat("\n============================================================\n")
cat("TOP-K SENSITIVITY ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("K values: 25, 50, 100\n")
cat("No RF/NPDR models were refit; saved full rankings were reused.\n")
cat("Bootstrap stability was NOT rerun and remains K=50.\n")
cat("Primary K=50 validation: PASSED\n")
cat("Main output:\n")
cat(file.path(tables_dir, "Supplementary_Table_S1_TopK_sensitivity_compact.csv"), "\n")
cat("============================================================\n\n")
print(compact, n = Inf)
