############################################################
## ARF_eval.R — End-to-end evaluation for ARF synthetic data
############################################################

## ========== 0) Paths & Setup ==========
set.seed(42)

# >>> EDIT THIS PATH IF NEEDED <<<
base_dir  <- "~/Desktop/My paper"
real_path <- file.path(base_dir, "gene_with_label.csv")
arf_path  <- file.path(base_dir, "synthetic_rna_arf.csv")
out_dir   <- file.path(base_dir, "ARF_eval")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========== 1) Packages ==========
need_pkgs <- c("tidyverse","ranger","umap","rsample","yardstick")
for(p in need_pkgs){
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it and re-run.", p))
  }
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(umap)
  library(rsample)
  library(yardstick)
})

venn_available <- requireNamespace("VennDiagram", quietly = TRUE)
npdr_available <- requireNamespace("npdr", quietly = TRUE)

# Helper to safely write tibbles
write_tbl <- function(x, path) readr::write_csv(tibble::as_tibble(x), path)

## ========== 2) Load data ==========
if(!file.exists(real_path)) stop("real_path not found: ", real_path)
if(!file.exists(arf_path))  stop("arf_path not found: ", arf_path)

real <- read.csv(real_path, check.names = FALSE)
arf  <- read.csv(arf_path,  check.names = FALSE)

if(!"class" %in% names(real)) stop("Real data must have a 'class' column.")
if(!"class" %in% names(arf))  stop("ARF data must have a 'class' column.")

# Map class to factor(HC/MDD) if numeric 0/1; otherwise coerce character to factor
if(is.numeric(real$class)) real$class <- factor(real$class, levels = c(0,1), labels = c("HC","MDD"))
if(is.numeric(arf$class))  arf$class  <- factor(arf$class,  levels = c(0,1), labels = c("HC","MDD"))
if(is.character(real$class) && !all(unique(real$class) %in% c("HC","MDD"))) real$class <- factor(real$class)
if(is.character(arf$class)  && !all(unique(arf$class)  %in% c("HC","MDD")))  arf$class  <- factor(arf$class)

## ========== 3) Harmonize & clean ==========
genes_real <- setdiff(names(real), "class")
genes_arf  <- setdiff(names(arf),  "class")
common     <- intersect(genes_real, genes_arf)

# Keep only shared genes + class
realX <- real[, c(common, "class")]
arfX  <- arf[,  c(common, "class")]

# Force numeric on gene columns (class remains factor)
realX[common] <- lapply(realX[common], \(x) as.numeric(as.character(x)))
arfX[common]  <- lapply(arfX[common],  \(x) as.numeric(as.character(x)))

# Drop zero-variance genes (based on Real)
nzv <- common[sapply(realX[common], sd, na.rm = TRUE) == 0]
if(length(nzv)){
  keep   <- setdiff(common, nzv)
  realX  <- realX[, c(keep,"class")]
  arfX   <- arfX[,  c(keep,"class")]
  common <- keep
}

# Remove rows with any NA
realX <- tidyr::drop_na(realX)
arfX  <- tidyr::drop_na(arfX)

## ========== 4) Utility: RF OOB accuracy (100 reps) ==========
rf_once <- function(df, ntree = 5000){
  m <- ranger(class ~ ., data = df, num.trees = ntree)
  as.numeric(1 - m$prediction.error)
}

R <- 100
acc_real <- replicate(R, rf_once(realX))
acc_arf  <- replicate(R, rf_once(arfX))
acc_comb <- replicate(R, {
  comb <- dplyr::bind_rows(
    dplyr::mutate(realX, .src = "Real"),
    dplyr::mutate(arfX,  .src = "ARF")
  ) |> dplyr::select(-.src)
  rf_once(comb)
})

acc_df <- tibble(replicate = 1:R, real = acc_real, arf = acc_arf, combined_arf = acc_comb)
write_tbl(acc_df, file.path(out_dir, "rf_accuracies_arf_vs_real_combined.csv"))

# Summary (mean ± sd) + paired t-test (Real+ARF vs Real)
acc_summary <- acc_df |>
  summarise(
    real_mean = mean(real), real_sd = sd(real),
    arf_mean  = mean(arf),  arf_sd  = sd(arf),
    comb_mean = mean(combined_arf), comb_sd = sd(combined_arf)
  )
print(acc_summary)
write_tbl(acc_summary, file.path(out_dir, "table_acc_summary_arf.csv"))

tt <- t.test(acc_df$combined_arf, acc_df$real, paired = TRUE)
dz <- as.numeric(tt$statistic) / sqrt(length(acc_df$real))  # Cohen's dz
eff_summary <- tibble(
  mean_diff = mean(acc_df$combined_arf - acc_df$real),
  ci_low    = tt$conf.int[1],
  ci_high   = tt$conf.int[2],
  t         = as.numeric(tt$statistic),
  df        = as.numeric(tt$parameter),
  p_value   = tt$p.value,
  cohens_dz = dz
)
print(eff_summary)
write_tbl(eff_summary, file.path(out_dir, "effectsize_real_vs_realplusarf.csv"))

## ========== 5) UMAP visualizations ==========
cfg <- umap.defaults; cfg$n_neighbors <- 15; cfg$min_dist <- 0.1

# (optional) for reproducible UMAP
set.seed(42)

# ARF-only
emb_arf <- as.data.frame(
  umap(as.matrix(arfX[, setdiff(names(arfX), "class")]), config = cfg)$layout
)
emb_arf$class <- arfX$class
p_arf <- ggplot(emb_arf, aes(V1, V2, color = class)) +
  geom_point(alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(plot.title = element_blank())
ggsave(file.path(out_dir, "umap_arf.png"), p_arf, width = 6, height = 5, dpi = 200)
p_arf
# Real vs ARF
comb_all <- dplyr::bind_rows(
  dplyr::mutate(realX, .src = "Real"),
  dplyr::mutate(arfX,  .src = "ARF")
)

emb_comb <- as.data.frame(
  umap(as.matrix(comb_all[, setdiff(names(comb_all), c("class", ".src"))]), config = cfg)$layout
)
emb_comb$src   <- comb_all$.src
emb_comb$class <- comb_all$class

p_comb <- ggplot(emb_comb, aes(V1, V2, color = src, shape = class)) +
  geom_point(alpha = 0.8) + theme_minimal() +
  labs(title = "UMAP - Real vs ARF")
ggsave(file.path(out_dir, "umap_real_vs_arf.png"), p_comb, width = 6, height = 5, dpi = 200)
p_comb

## ========== 6) Fidelity: KS + Wilcoxon + FDR ==========
ks_p <- sapply(common, \(g){
  x <- realX[[g]]; y <- arfX[[g]]
  if(all(is.finite(x)) && all(is.finite(y))) suppressWarnings(ks.test(x,y)$p.value) else NA_real_
})
wil_p <- sapply(common, \(g){
  x <- realX[[g]]; y <- arfX[[g]]
  if(all(is.finite(x)) && all(is.finite(y))) suppressWarnings(wilcox.test(x,y)$p.value) else NA_real_
})

dist_df <- tibble(gene = common,
                  ks_pvalue = ks_p,
                  wil_pvalue = wil_p) |>
  mutate(q_ks  = p.adjust(ks_pvalue,  method = "BH"),
         q_wil = p.adjust(wil_pvalue, method = "BH"))
write_tbl(dist_df, file.path(out_dir, "distribution_checks_real_vs_arf_FDR.csv"))

fidelity_summary <- dist_df |>
  summarise(
    ks_sig_pct_raw  = mean(ks_pvalue  < 0.05, na.rm = TRUE) * 100,
    wil_sig_pct_raw = mean(wil_pvalue < 0.05, na.rm = TRUE) * 100,
    ks_FDR_sig_pct  = mean(q_ks       < 0.05, na.rm = TRUE) * 100,
    wil_FDR_sig_pct = mean(q_wil      < 0.05, na.rm = TRUE) * 100
  )
print(fidelity_summary)
write_tbl(fidelity_summary, file.path(out_dir, "fidelity_summary_arf.csv"))

## ========== 7) Privacy: domain-separator AUC (5×5 CV) ==========
# Goal: can a classifier tell "Real" vs "ARF"? Closer to 0.5 → better privacy
comb_bin <- comb_all
comb_bin$src_bin <- factor(ifelse(comb_bin$.src == "Real", "Real", "ARF"))

folds <- vfold_cv(comb_bin, v = 5, repeats = 5, strata = src_bin)
oof <- lapply(folds$splits, function(s){
  train <- analysis(s); test <- assessment(s)
  m <- ranger(src_bin ~ . - .src - class, data = train, probability = TRUE, num.trees = 2000)
  pr <- predict(m, test)$predictions[, "Real"]
  tibble(truth = test$src_bin, .pred_Real = pr)
}) |> bind_rows()

auc_cv <- yardstick::roc_auc(oof, truth = truth, .pred_Real)
print(auc_cv)
write_tbl(as.data.frame(auc_cv), file.path(out_dir, "privacy_auc_cv_arf.csv"))

# Also report "adjusted-to-≥0.5" convention (so 0.40 -> 0.60)
auc_val    <- as.numeric(auc_cv$.estimate)
auc_report <- ifelse(auc_val < 0.5, 1 - auc_val, auc_val)
write_tbl(tibble(auc_raw = auc_val, auc_report = auc_report),
          file.path(out_dir, "privacy_auc_cv_arf_normalized.csv"))

## ========== 8) Feature importance — RF (permutation) ==========
TopK <- 50  # change if needed

rf_imp_tbl <- function(df){
  m <- ranger(
    class ~ ., data = df,
    num.trees = 5000,
    importance = "permutation",
    scale.permutation.importance = TRUE
  )
  tibble(gene = names(m$variable.importance),
         rf_importance = as.numeric(m$variable.importance)) |>
    arrange(desc(rf_importance))
}

imp_real_tbl <- rf_imp_tbl(realX)
imp_arf_tbl  <- rf_imp_tbl(arfX)
write_tbl(imp_real_tbl, file.path(out_dir, "rf_importance_real.csv"))
write_tbl(imp_arf_tbl,  file.path(out_dir, "rf_importance_arf.csv"))

# Top-K + overlap + Jaccard
top_real_rf <- head(imp_real_tbl$gene, TopK)
top_arf_rf  <- head(imp_arf_tbl$gene,  TopK)
overlap_rf  <- intersect(top_real_rf, top_arf_rf)
jaccard_rf  <- length(overlap_rf) / length(union(top_real_rf, top_arf_rf))

rf_overlap_sum <- tibble(
  TopK      = TopK,
  n_real    = length(top_real_rf),
  n_arf     = length(top_arf_rf),
  n_overlap = length(overlap_rf),
  jaccard   = jaccard_rf
)
print(rf_overlap_sum)
write_tbl(rf_overlap_sum, file.path(out_dir, "feature_overlap_rf_summary.csv"))
write_tbl(tibble(gene = overlap_rf), file.path(out_dir, "feature_overlap_rf_genes.csv"))

# Optional barplots
TopBar <- 20
p_real_bar <- imp_real_tbl |> slice_head(n = min(TopBar, nrow(imp_real_tbl))) |>
  ggplot(aes(x = reorder(gene, rf_importance), y = rf_importance)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(title = sprintf("RF Importance (Real) - Top %d", TopBar),
       x = "Gene", y = "Permutation Importance")
ggsave(file.path(out_dir, "rf_importance_real_top20.png"), p_real_bar, width = 7, height = 6, dpi = 200)

p_arf_bar <- imp_arf_tbl |> slice_head(n = min(TopBar, nrow(imp_arf_tbl))) |>
  ggplot(aes(x = reorder(gene, rf_importance), y = rf_importance)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(title = sprintf("RF Importance (ARF) - Top %d", TopBar),
       x = "Gene", y = "Permutation Importance")
ggsave(file.path(out_dir, "rf_importance_arf_top20.png"), p_arf_bar, width = 7, height = 6, dpi = 200)

# Optional Venn (RF)
if (venn_available) {
  suppressPackageStartupMessages({ library(VennDiagram); library(grid) })
  venn_sets <- list(Real = top_real_rf); venn_sets[["ARF"]] <- top_arf_rf
  venn_plot <- VennDiagram::venn.diagram(
    x = venn_sets,
    filename = NULL, output = TRUE,
    fill = c("#8da0cb", "#fc8d62"), alpha = 0.6,
    cex = 1.2, cat.cex = 1.2, lwd = 1.5
  )
  png(file.path(out_dir, sprintf("venn_rf_top%d_real_vs_arf.png", TopK)),
      width = 900, height = 650, res = 120)
  grid::grid.draw(venn_plot); dev.off()
}

## ========== 9) Feature importance — NPDR (Real / ARF / Real+ARF) ==========
if (npdr_available) {
  suppressPackageStartupMessages(library(npdr))
  
  # ---- Helpers (prefer adjusted p-value; fallback to |beta|) ----
  pick_num_vec <- function(df, candidates){
    nm <- intersect(candidates, names(df))
    if (length(nm) == 0) return(rep(NA_real_, nrow(df)))
    v  <- suppressWarnings(as.numeric(df[[nm[1]]]))
    if (length(v) == 0) v <- rep(NA_real_, nrow(df))
    v
  }
  rank_npdr <- function(tbl){
    p_adj <- pick_num_vec(tbl, c("pval.adj","padj","p_adj","adj.p","adj.P.Val","q.value","qval","q_adj"))
    beta  <- pick_num_vec(tbl, c("beta","coef","coefficient","importance","score","effect.size","effect"))
    p_use <- ifelse(is.finite(p_adj), p_adj, Inf)
    b_use <- ifelse(is.finite(beta),  beta,  0)
    ord   <- order(p_use, -abs(b_use), seq_len(nrow(tbl)))
    out   <- tibble::as_tibble(tbl)[ord, , drop = FALSE]
    out$p_adj_use <- p_use[ord]; out$beta_use <- b_use[ord]
    out
  }
  run_npdr_tbl <- function(df, label){
    X <- scale(as.matrix(dplyr::select(df, -class)))
    y <- droplevels(df$class)
    mod <- npdr::npdr(
      outcome         = y,
      dataset         = X,
      regression.type = "binomial",
      nbd.method      = "multisurf",
      nbd.metric      = "euclidean",
      knn             = 20,
      padj.method     = "fdr",
      verbose         = FALSE
    )
    tbl <- tibble::as_tibble(mod)
    if (!"gene" %in% names(tbl)) {
      if (nrow(tbl) == ncol(X)) tbl$gene <- colnames(X) else tbl$gene <- paste0("gene_", seq_len(nrow(tbl)))
    }
    ranked <- rank_npdr(tbl)
    readr::write_csv(ranked, file.path(out_dir, paste0("npdr_importance_", tolower(label), ".csv")))
    ranked
  }
  
  # --- Run NPDR on three sets ---
  npdr_real_rank <- run_npdr_tbl(realX, "Real")
  npdr_arf_rank  <- run_npdr_tbl(arfX,  "ARF")
  npdr_comb_rank <- run_npdr_tbl(dplyr::bind_rows(realX, arfX), "RealPlusARF")
  
  # --- Pick TopK (prefer q<0.05; else TopK with independent tie-breaks) ---
  select_top <- function(ranked_tbl, TopK){
    sig <- dplyr::filter(ranked_tbl, is.finite(p_adj_use) & p_adj_use < 0.05)
    if (nrow(sig) > 0) {
      list(genes = head(sig[order(sig$p_adj_use, -abs(sig$beta_use)), "gene", drop=TRUE], min(TopK, nrow(sig))),
           note  = "significant (q<0.05)")
    } else {
      set.seed(1000 + nrow(ranked_tbl))
      tmp <- ranked_tbl |>
        dplyr::arrange(p_adj_use, dplyr::desc(abs(beta_use))) |>
        dplyr::mutate(tie = runif(dplyr::n())) |>
        dplyr::arrange(p_adj_use, dplyr::desc(abs(beta_use)), tie)
      list(genes = head(tmp$gene, TopK),
           note  = "no significant genes; TopK with random tie-break")
    }
  }
  
  sel_real <- select_top(npdr_real_rank, TopK)
  sel_arf  <- select_top(npdr_arf_rank,  TopK)
  sel_comb <- select_top(npdr_comb_rank, TopK)
  
  top_real_npdr <- sel_real$genes
  top_arf_npdr  <- sel_arf$genes
  top_comb_npdr <- sel_comb$genes
  
  # --- Pairwise overlaps & Jaccard ---
  jac <- function(a,b) length(intersect(a,b)) / length(union(a,b))
  pair_summ <- tibble::tibble(
    pair     = c("Real vs ARF", "Real vs Real+ARF", "ARF vs Real+ARF"),
    TopK_eff = c(length(top_real_npdr), length(top_real_npdr), length(top_arf_npdr)),
    n_A      = c(length(top_real_npdr), length(top_real_npdr), length(top_arf_npdr)),
    n_B      = c(length(top_arf_npdr),  length(top_comb_npdr), length(top_comb_npdr)),
    n_overlap= c(length(intersect(top_real_npdr, top_arf_npdr)),
                 length(intersect(top_real_npdr, top_comb_npdr)),
                 length(intersect(top_arf_npdr,  top_comb_npdr))),
    jaccard  = c(jac(top_real_npdr, top_arf_npdr),
                 jac(top_real_npdr, top_comb_npdr),
                 jac(top_arf_npdr,  top_comb_npdr))
  )
  print(pair_summ)
  readr::write_csv(pair_summ, file.path(out_dir, "npdr_pairwise_overlap_summary_real_arf_realplusarf.csv"))
  
  # Optional Venns
  if (venn_available) {
    suppressPackageStartupMessages({library(VennDiagram); library(grid)})
    # Real vs Real+ARF
    v1_sets <- list(Real = top_real_npdr, `Real+ARF` = top_comb_npdr)
    v1 <- VennDiagram::venn.diagram(x = v1_sets, filename = NULL, output = TRUE,
                                    fill = c("#8da0cb","#66c2a5"), alpha = 0.6,
                                    cex = 1.2, cat.cex = 1.2, lwd = 1.5)
    png(file.path(out_dir, sprintf("venn_npdr_top%d_real_vs_realplusarf.png", length(top_real_npdr))),
        width = 900, height = 650, res = 120); grid::grid.draw(v1); dev.off()
    
    # 3-way (Real, ARF, Real+ARF)
    v3_sets <- list(Real = top_real_npdr, ARF = top_arf_npdr, `Real+ARF` = top_comb_npdr)
    v3 <- VennDiagram::venn.diagram(x = v3_sets, filename = NULL, output = TRUE,
                                    fill = c("#8da0cb","#fc8d62","#66c2a5"), alpha = 0.6,
                                    cex = 1.0, cat.cex = 1.1, lwd = 1.5)
    png(file.path(out_dir, sprintf("venn_npdr_top%d_real_arf_realplusarf.png", TopK)),
        width = 950, height = 700, res = 120); grid::grid.draw(v3); dev.off()
  }
} else {
  message("Package 'npdr' not installed; skipping NPDR section.")
}

## ========== 10) Dilution metric: Real Top-K rank drop after adding ARF ==========
TopK <- min(50, nrow(realX)-1)  # safety

rf_imp_tbl2 <- function(df){
  m <- ranger(class ~ ., data = df, num.trees = 3000,
              importance = "permutation", scale.permutation.importance = TRUE)
  tibble(gene = names(m$variable.importance),
         imp  = as.numeric(m$variable.importance)) |>
    arrange(desc(imp))
}
imp_real    <- rf_imp_tbl2(realX)
imp_realplA <- rf_imp_tbl2(dplyr::bind_rows(realX, arfX))

# Ranks (*** FIX: rename before join to avoid 'rank_realplus' missing ***)
rank_tbl <- function(tbl) tbl |> mutate(rank = rank(-imp, ties.method = "average")) |> select(gene, rank)
r_real    <- rank_tbl(imp_real)    |> rename(rank_real     = rank)
r_realplA <- rank_tbl(imp_realplA) |> rename(rank_realplus = rank)

top_real <- head(imp_real$gene, TopK)
dilute_arf <- r_real |>
  filter(gene %in% top_real) |>
  left_join(r_realplA, by = "gene") |>
  mutate(
    delta_rank = rank_realplus - rank_real,
    dropped    = as.integer(rank_realplus > TopK)
  )

dilute_summary_arf <- dilute_arf |>
  summarise(
    mean_delta   = mean(delta_rank, na.rm = TRUE),
    median_delta = median(delta_rank, na.rm = TRUE),
    pct_dropped  = mean(dropped, na.rm = TRUE) * 100
  )
print(dilute_summary_arf)
readr::write_csv(dilute_arf,         file.path(out_dir, "dilution_detail_real_to_realplusARF_RF.csv"))
readr::write_csv(dilute_summary_arf, file.path(out_dir, "dilution_summary_real_to_realplusARF_RF.csv"))

## ========== 11) Stability of top-K with bootstraps (Jaccard distributions) ==========
set.seed(123)
B <- 50
TopK <- 50

get_topk_rf <- function(df, k=TopK) {
  m <- ranger(class ~ ., data=df, num.trees=2000,
              importance="permutation", scale.permutation.importance=TRUE)
  tibble(gene = names(m$variable.importance),
         imp  = as.numeric(m$variable.importance)) |>
    arrange(desc(imp)) |>
    slice_head(n = k) |>
    pull(gene)
}
jac <- function(a,b) length(intersect(a,b))/length(union(a,b))

boot_jaccards <- function(dfA, dfB, B=50, k=TopK) {
  out <- vector("numeric", B)
  nA  <- nrow(dfA); nB <- nrow(dfB)
  for (b in seq_len(B)) {
    A_bs <- dfA[sample.int(nA, replace = TRUE), ]
    B_bs <- dfB[sample.int(nB, replace = TRUE), ]
    topA <- get_topk_rf(A_bs, k)
    topB <- get_topk_rf(B_bs, k)
    out[b] <- jac(topA, topB)
  }
  out
}

# Cross-set stability: Real vs ARF
jac_real_arf <- boot_jaccards(realX, arfX, B=B, k=TopK)

# Within-set self-consistency
boot_self <- function(df, B=50, k=TopK) {
  out <- vector("numeric", B); n <- nrow(df)
  for (b in seq_len(B)) {
    A <- df[sample.int(n, replace = TRUE), ]
    B2<- df[sample.int(n, replace = TRUE), ]
    out[b] <- jac(get_topk_rf(A,k), get_topk_rf(B2,k))
  }
  out
}
stab_real <- boot_self(realX, B=B, k=TopK)
stab_arf  <- boot_self(arfX,  B=B, k=TopK)

tibble(
  set = c(rep("Real_vs_ARF", length(jac_real_arf)),
          rep("Real_self",   length(stab_real)),
          rep("ARF_self",    length(stab_arf))),
  jaccard = c(jac_real_arf, stab_real, stab_arf)
) |>
  readr::write_csv(file.path(out_dir, "jaccard_bootstrap_RF.csv"))

############################################################
## End of ARF_eval.R
############################################################
