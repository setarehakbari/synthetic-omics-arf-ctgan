############################################################
# make-figures.R  (clean UMAP + all figures & tables)
# Output:  ~/Desktop/My paper/results
############################################################

## ========== 0) Setup ==========
base_dir <- path.expand("~/Desktop/My paper")
dir_arf  <- file.path(base_dir, "ARF_eval")
dir_ctg  <- file.path(base_dir, "CTGAN_eval")
res_dir  <- file.path(base_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# (اختیاری) مسیرِ داده‌های خام برای ساخت UMAP تمیز
real_csv <- file.path(base_dir, "gene_with_label.csv")
arf_csv  <- file.path(base_dir, "synthetic_rna_arf.csv")
ctg_csv  <- file.path(base_dir, "synthetic_rna_ctgan.csv")

set.seed(42)

## ========== 1) Packages ==========
need <- c("tidyverse","cowplot","scales","glue","patchwork")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) stop(sprintf("Install '%s'", p))
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(scales); library(glue); library(patchwork)
})

has_umap   <- requireNamespace("umap", quietly = TRUE)
has_gt     <- requireNamespace("gt", quietly = TRUE)
has_webshot<- requireNamespace("webshot2", quietly = TRUE)

## ========== Helpers ==========
read_tbl <- function(path){
  if (file.exists(path)) suppressWarnings(readr::read_csv(path, show_col_types = FALSE)) else { message("Missing: ", path); NULL }
}
save_figure <- function(plot_obj, name, width=8, height=5, dpi=300){
  ggsave(file.path(res_dir, paste0(name,".png")), plot_obj, width=width, height=height, dpi=dpi)
  # PDF با cairo اگر موجود بود
  pdf_path <- file.path(res_dir, paste0(name,".pdf"))
  ok <- TRUE
  tryCatch({ ggsave(pdf_path, plot_obj, width=width, height=height, device=cairo_pdf) },
           error=function(e){ ok <<- FALSE })
  if(!ok) ggsave(pdf_path, plot_obj, width=width, height=height, device="pdf")
}
theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        legend.position = "right")

## =========================================================
## FIGURE 1 — Utility (RF OOB accuracy)
## =========================================================
acc_arf <- read_tbl(file.path(dir_arf, "rf_accuracies_arf_vs_real_combined.csv")) |>
  mutate(model="ARF") |> rename(synth=arf, combined=combined_arf)
acc_ctg <- read_tbl(file.path(dir_ctg,"rf_accuracies_ctgan_vs_real_combined.csv")) |>
  mutate(model="CTGAN") |> rename(synth=ctgan, combined=combined_ctgan)

if (!is.null(acc_arf) && !is.null(acc_ctg)) {
  acc_long <- bind_rows(acc_arf, acc_ctg) |>
    pivot_longer(c(real,synth,combined), names_to="set", values_to="accuracy") |>
    mutate(set=factor(set, levels=c("real","synth","combined"),
                      labels=c("Real","Synthetic","Real+Synthetic")))
  p_acc <- ggplot(acc_long, aes(set, accuracy)) +
    geom_violin(trim=FALSE, alpha=.9) +
    geom_boxplot(width=.16, outlier.shape=NA) +
    stat_summary(fun=mean, geom="point", shape=4, size=3, stroke=1.1) +
    facet_wrap(~model, nrow=1) +
    coord_cartesian(ylim=c(0.45,0.70)) +
    labs(title="Utility: RF OOB accuracy (100 reps)", x=NULL, y="Accuracy (OOB)") +
    theme_pub
  save_figure(p_acc, "Fig1_utility_boxplots", 10, 5)
}

## =========================================================
## FIGURE 2 — UMAP (clean build from CSV if available; else fallback montage)
## =========================================================
## =========================================================
## FIGURE 2 — UMAP grids (clean titles, no overlap)
## =========================================================
png_paths <- list(
  ARF   = c(file.path(dir_arf,"umap_real_vs_arf.png"),
            file.path(dir_arf,"umap_arf.png")),
  CTGAN = c(file.path(dir_ctg,"umap_real_vs_ctgan.png"),
            file.path(dir_ctg,"umap_ctgan.png"))
)

# build a single tile with image + subtitle placed OUTSIDE the panel
build_umap_tile <- function(path, subtitle) {
  if (!file.exists(path)) {
    return(ggdraw() + draw_label(glue::glue("Missing: {basename(path)}"), size=12))
  }
  g <- ggdraw() + 
    draw_image(path, scale = 1) +
    theme(plot.margin = margin(t = 8, r = 8, b = 8, l = 8))
  cowplot::add_sub(g, subtitle, x = 0, hjust = 0, vpadding = grid::unit(6, "pt"))
}

row_arf <- plot_grid(
  build_umap_tile(png_paths$ARF[1], "Real vs ARF (UMAP)"),
  build_umap_tile(png_paths$ARF[2], "ARF-only (UMAP)"),
  nrow = 1, rel_widths = c(1,1)
)

row_ctg <- plot_grid(
  build_umap_tile(png_paths$CTGAN[1], "Real vs CTGAN (UMAP)"),
  build_umap_tile(png_paths$CTGAN[2], "CTGAN-only (UMAP)"),
  nrow = 1, rel_widths = c(1,1)
)

# add a small spacer between rows so subtitles don't touch
p_umap <- plot_grid(
  ggdraw() + draw_label("UMAP embeddings", fontface="bold", size=14),
  row_arf,
  ggdraw(),  # spacer
  row_ctg,
  ncol = 1,
  rel_heights = c(0.10, 1, 0.06, 1)
)
save_figure(p_umap, "Fig2_umap_grids", width = 11, height = 9)  # کمی پهن‌تر برای جا داشتن زیرنویس‌ها


## =========================================================
## FIGURE 3 — Privacy AUC
## =========================================================
auc_arf_raw <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf.csv"))
auc_arf_adj <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf_normalized.csv"))
auc_ctg_raw <- read_tbl(file.path(dir_ctg,"privacy_auc_cv_ctgan.csv"))
auc_ctg_adj <- read_tbl(file.path(dir_ctg,"privacy_auc_cv_ctgan_normalized.csv"))

if (!is.null(auc_arf_raw) & !is.null(auc_ctg_raw) & !is.null(auc_arf_adj) & !is.null(auc_ctg_adj)) {
  df_priv <- tibble(
    model=c("ARF","CTGAN"),
    AUC_raw = c(auc_arf_raw$.estimate[1], auc_ctg_raw$.estimate[1]),
    AUC_norm= c(auc_arf_adj$auc_report[1], auc_ctg_adj$auc_report[1])
  ) |> pivot_longer(c(AUC_raw,AUC_norm), names_to="type", values_to="AUC")
  p_priv <- ggplot(df_priv, aes(model, AUC, fill=type)) +
    geom_col(position=position_dodge(.6), width=.55) +
    geom_hline(yintercept=.5, linetype=2) +
    scale_y_continuous(limits=c(0,1)) +
    scale_fill_manual(values=c("#9ecae1","#6baed6"), labels=c("Raw AUC","Normalized (≥0.5)")) +
    labs(title="Privacy: Domain-separator AUC (5×5 CV)", x=NULL, y="AUC") +
    theme_pub + theme(legend.position="top")
  save_figure(p_priv, "Fig3_privacy_auc", 7, 5)
}

## =========================================================
## FIGURE 4 — Fidelity summary
## =========================================================
fid_arf <- read_tbl(file.path(dir_arf,"fidelity_summary_arf.csv"))
fid_ctg <- read_tbl(file.path(dir_ctg,"fidelity_summary_ctgan.csv"))
if (!is.null(fid_arf) & !is.null(fid_ctg)) {
  df_fid <- bind_rows(fid_arf |> mutate(model="ARF"),
                      fid_ctg |> mutate(model="CTGAN")) |>
    select(model, ks_FDR_sig_pct, wil_FDR_sig_pct) |>
    pivot_longer(-model, names_to="test", values_to="pct") |>
    mutate(test=recode(test, ks_FDR_sig_pct="KS (BH-FDR)", wil_FDR_sig_pct="Wilcoxon (BH-FDR)"))
  p_fid <- ggplot(df_fid, aes(test, pct, fill=model)) +
    geom_col(position=position_dodge(.6), width=.55) +
    scale_y_continuous(labels=percent_format(scale=1), limits=c(0,100)) +
    labs(title="Fidelity: % genes significant (BH-FDR)", x=NULL, y="% significant") +
    theme_pub + theme(legend.position="top")
  save_figure(p_fid, "Fig4_fidelity_summary", 8, 5)
}

## =========================================================
## FIGURE 5 — RF Top-20
## =========================================================
rf_real_arf <- read_tbl(file.path(dir_arf,"rf_importance_real.csv"))
rf_arf      <- read_tbl(file.path(dir_arf,"rf_importance_arf.csv"))
rf_real_ctg <- read_tbl(file.path(dir_ctg,"rf_importance_real.csv"))
rf_ctg      <- read_tbl(file.path(dir_ctg,"rf_importance_ctgan.csv"))

make_top20 <- function(tbl, title, value_col="rf_importance"){
  if (is.null(tbl)) return(ggdraw() + draw_label(glue("Missing: {title}")))
  tbl |> slice_head(n=min(20,nrow(tbl))) |>
    mutate(gene=fct_reorder(gene, !!sym(value_col))) |>
    ggplot(aes(gene, !!sym(value_col))) + geom_col() + coord_flip() +
    labs(title=title, x="Gene", y="Permutation Importance") + theme_pub
}
p_rf <- plot_grid(
  ggdraw() + draw_label("Top-20 RF permutation importance", fontface="bold", size=14),
  plot_grid(make_top20(rf_real_arf,"RF Importance (Real) — ARF run"),
            make_top20(rf_arf,"RF Importance (ARF)"), nrow=1),
  plot_grid(make_top20(rf_real_ctg,"RF Importance (Real) — CTGAN run"),
            make_top20(rf_ctg,"RF Importance (CTGAN)"), nrow=1),
  ncol=1, rel_heights=c(.12,1,1)
)
save_figure(p_rf, "Fig5_rf_importance_top20", 12, 10)

## =========================================================
## FIGURE 6 — NPDR overlaps
## =========================================================
npdr_arf <- read_tbl(file.path(dir_arf,"npdr_pairwise_overlap_summary_real_arf_realplusarf.csv"))
npdr_ctg <- read_tbl(file.path(dir_ctg,"npdr_pairwise_overlap_summary_real_ctgan_realplusctgan.csv"))

plot_npdr_bar <- function(df, title){
  if (is.null(df)) return(ggdraw() + draw_label(glue("Missing: {title}")))
  df |> mutate(pair=factor(pair, levels=pair)) |>
    ggplot(aes(pair, jaccard)) +
    geom_col(width=.6) +
    geom_text(aes(label=sprintf("n=%d", n_overlap)), vjust=-.35, size=3.4) +
    scale_y_continuous(limits=c(0,1)) +
    labs(title=title, x=NULL, y="Jaccard") + theme_pub
}
p_npdr <- plot_grid(plot_npdr_bar(npdr_arf,"NPDR Top-K overlap — ARF"),
                    plot_npdr_bar(npdr_ctg,"NPDR Top-K overlap — CTGAN"), nrow=1)
save_figure(p_npdr, "Fig6_npdr_overlap_jaccard", 10, 4.5)

## =========================================================
## FIGURE 7 — Dilution
## =========================================================
dil_arf <- read_tbl(file.path(dir_arf,"dilution_detail_real_to_realplusARF_RF.csv")) |> mutate(model="ARF")
dil_ctg <- read_tbl(file.path(dir_ctg,"dilution_detail_real_to_realplusCTGAN_RF.csv")) |> mutate(model="CTGAN")
if (!is.null(dil_arf) & !is.null(dil_ctg)) {
  dil_all <- bind_rows(dil_arf, dil_ctg)
  p_d1 <- ggplot(dil_all, aes(delta_rank, fill=model)) +
    geom_histogram(bins=40, alpha=.9, position="identity") +
    labs(title="Dilution: Δ rank of Real Top-K after adding Synthetic", x="Δ rank", y="Count") +
    theme_pub
  p_d2 <- dil_all |> group_by(model) |>
    summarise(mean_delta=mean(delta_rank,na.rm=TRUE),
              median_delta=median(delta_rank,na.rm=TRUE),
              pct_dropped=mean(dropped==1,na.rm=TRUE)*100) |>
    pivot_longer(-model, names_to="metric", values_to="value") |>
    ggplot(aes(metric, value, fill=model)) +
    geom_col(position=position_dodge(.7), width=.6) +
    labs(title="Dilution summary", x=NULL, y="Value") +
    theme_pub + theme(axis.text.x=element_text(angle=18, hjust=1))
  save_figure(plot_grid(p_d1,p_d2,ncol=1,rel_heights=c(1,.9)), "Fig7_dilution_metric", 10, 9)
}

## =========================================================
## FIGURE 8 — Stability (bootstrap Jaccard)
## =========================================================
stab_arf <- read_tbl(file.path(dir_arf,"jaccard_bootstrap_RF.csv")) |> mutate(model="ARF")
stab_ctg <- read_tbl(file.path(dir_ctg,"jaccard_bootstrap_RF.csv")) |> mutate(model="CTGAN")
if (!is.null(stab_arf) & !is.null(stab_ctg)) {
  stab_all <- bind_rows(stab_arf, stab_ctg) |>
    mutate(set=factor(set, levels=c("Real_self","ARF_self","Real_vs_ARF","CTGAN_self","Real_vs_CTGAN")))
  p_stab <- ggplot(stab_all, aes(set, jaccard, fill=model)) +
    geom_violin(trim=FALSE, alpha=.9) +
    geom_boxplot(width=.16, outlier.shape=NA, position=position_dodge(.75)) +
    coord_cartesian(ylim=c(0,1)) +
    labs(title="Stability of Top-K (Jaccard over bootstraps)", x=NULL, y="Jaccard") +
    theme_pub + theme(axis.text.x=element_text(angle=18,hjust=1))
  save_figure(p_stab, "Fig8_stability_jaccard", 10, 5.5)
}

## ========================= Tables (1–5) =========================
save_pretty_table <- function(df, title, stem, digits=4){
  out_csv <- file.path(res_dir, paste0(stem,".csv")); readr::write_csv(df, out_csv)
  if (has_gt && has_webshot) {
    library(gt)
    tb <- df |> gt() |> tab_header(title=md(paste0("**",title,"**"))) |>
      fmt_number(columns=where(is.numeric), decimals=digits)
    gt::gtsave(tb, file.path(res_dir, paste0(stem,".png")))
    gt::gtsave(tb, file.path(res_dir, paste0(stem,".pdf")))
  }
}

# Table 1 — Utility summary + paired t-test
arf_sum <- read_tbl(file.path(dir_arf,"table_acc_summary_arf.csv")) |> mutate(model="ARF")
arf_t   <- read_tbl(file.path(dir_arf,"effectsize_real_vs_realplusarf.csv")) |> mutate(model="ARF")
ctg_sum <- read_tbl(file.path(dir_ctg,"table_acc_summary_ctgan.csv")) |> mutate(model="CTGAN")
ctg_t   <- read_tbl(file.path(dir_ctg,"effectsize_real_vs_realplusctgan.csv")) |> mutate(model="CTGAN")
if (!is.null(arf_sum) & !is.null(arf_t) & !is.null(ctg_sum) & !is.null(ctg_t)) {
  norm_cols <- function(df){
    if ("arf_mean" %in% names(df)) df |> rename(synth_mean=arf_mean, synth_sd=arf_sd)
    else if ("ctg_mean" %in% names(df)) df |> rename(synth_mean=ctg_mean, synth_sd=ctg_sd)
    else df
  }
  tab1 <- bind_rows(left_join(norm_cols(arf_sum), arf_t, by="model"),
                    left_join(norm_cols(ctg_sum), ctg_t, by="model")) |>
    select(model, real_mean, real_sd, synth_mean, synth_sd, comb_mean, comb_sd,
           mean_diff, ci_low, ci_high, t, df, p_value, cohens_dz)
  save_pretty_table(tab1, "Table 1. Utility summary and paired t-tests", "Table1_utility_summary_with_ttest", 4)
}

# Table 2 — Privacy AUC
if (!is.null(auc_arf_raw) & !is.null(auc_ctg_raw) & !is.null(auc_arf_adj) & !is.null(auc_ctg_adj)) {
  tab2 <- tibble(model=c("ARF","CTGAN"),
                 AUC_raw=c(auc_arf_raw$.estimate[1], auc_ctg_raw$.estimate[1]),
                 AUC_norm=c(auc_arf_adj$auc_report[1], auc_ctg_adj$auc_report[1]))
  save_pretty_table(tab2, "Table 2. Privacy via domain-separator AUC (5×5 CV)", "Table2_privacy_auc", 3)
}

# Table 3 — Fidelity summary
if (!is.null(fid_arf) & !is.null(fid_ctg)) {
  tab3 <- bind_rows(fid_arf |> mutate(model="ARF"),
                    fid_ctg |> mutate(model="CTGAN")) |>
    select(model, ks_sig_pct_raw, wil_sig_pct_raw, ks_FDR_sig_pct, wil_FDR_sig_pct)
  save_pretty_table(tab3, "Table 3. Fidelity summary (% significant; raw & BH-FDR)", "Table3_fidelity_summary", 1)
}

# Table 4 — RF overlap
rf_ov_arf <- read_tbl(file.path(dir_arf,"feature_overlap_rf_summary.csv")) |> mutate(model="ARF")
rf_ov_ctg <- read_tbl(file.path(dir_ctg,"feature_overlap_rf_summary.csv")) |> mutate(model="CTGAN")
if (!is.null(rf_ov_arf) & !is.null(rf_ov_ctg)) {
  tab4 <- bind_rows(rf_ov_arf, rf_ov_ctg) |>
    select(model, TopK, n_real, n_overlap, jaccard,
           n_synth = any_of(c("n_arf","n_ctg")))
  save_pretty_table(tab4, "Table 4. RF Top-K overlap & Jaccard (Real vs Synthetic)", "Table4_rf_overlap", 3)
}

# Table 5 — Dilution summary
dil_sum_arf <- read_tbl(file.path(dir_arf,"dilution_summary_real_to_realplusARF_RF.csv")) |> mutate(model="ARF")
dil_sum_ctg <- read_tbl(file.path(dir_ctg,"dilution_summary_real_to_realplusCTGAN_RF.csv")) |> mutate(model="CTGAN")
if (!is.null(dil_sum_arf) & !is.null(dil_sum_ctg)) {
  tab5 <- bind_rows(dil_sum_arf, dil_sum_ctg) |> select(model, mean_delta, median_delta, pct_dropped)
  save_pretty_table(tab5, "Table 5. Dilution summary of Real Top-K after adding Synthetic", "Table5_dilution_summary", 1)
}

message("All done. Results saved to: ", res_dir)
############################################################
# End
############################################################
