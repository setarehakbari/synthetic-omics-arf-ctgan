############################################################
# make-figures.R
# Build publication-ready figures from outputs of:
############################################################

## ========== 0) Setup ==========
base_dir <- path.expand("~/Desktop/My paper")
dir_arf  <- file.path(base_dir, "ARF_eval")
dir_ctg  <- file.path(base_dir, "CTGAN_eval")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ========== 1) Packages ==========
need <- c("tidyverse","cowplot","scales","glue","patchwork")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it and re-run.", p))
  }
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(scales)
  library(glue)
  library(patchwork)
})

## Helpers
read_tbl <- function(path) {
  if (file.exists(path)) {
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  } else {
    message("Missing file, skipping: ", path)
    NULL
  }
}
save_figure <- function(plot_obj, name, width=8, height=5, dpi=300) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot_obj, width=width, height=height, dpi=dpi)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot_obj, width=width, height=height, device=cairo_pdf)
}

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face="bold"),
    strip.text       = element_text(face="bold")
  )

## ========== 2) FIGURE 1 — Utility (RF OOB accuracy distributions) ==========
# ARF
acc_arf <- read_tbl(file.path(dir_arf, "rf_accuracies_arf_vs_real_combined.csv")) |>
  mutate(model = "ARF") |> rename(synth = arf, combined = combined_arf)

# CTGAN
acc_ctg <- read_tbl(file.path(dir_ctg, "rf_accuracies_ctgan_vs_real_combined.csv")) |>
  mutate(model = "CTGAN") |> rename(synth = ctgan, combined = combined_ctgan)

if (!is.null(acc_arf) && !is.null(acc_ctg)) {
  acc_long <- bind_rows(acc_arf, acc_ctg) |>
    pivot_longer(cols = c(real, synth, combined),
                 names_to = "set", values_to = "accuracy") |>
    mutate(
      set = factor(set, levels = c("real","synth","combined"),
                   labels = c("Real","Synthetic","Real+Synthetic")),
      model = factor(model, levels = c("ARF","CTGAN"))
    )
  
  dodge <- position_dodge(width = 0.75)
  
  p_acc <- ggplot(acc_long, aes(x = set, y = accuracy, fill = model)) +
    geom_violin(
      trim = FALSE, alpha = 0.75,
      position = dodge
    ) +
    geom_boxplot(
      width = 0.18, outlier.shape = NA,
      position = dodge
    ) +
    stat_summary(
      fun = mean, geom = "point",
      position = dodge, shape = 4, size = 3, stroke = 1.1
    ) +
    coord_cartesian(ylim = c(0.45, 0.70)) +
    labs(
      title = "Utility: RF OOB accuracy (100 reps)",
      x = NULL, y = "Accuracy (OOB)", fill = "Generator"
    ) +
    theme_pub +
    theme(legend.position = "top")
  
  save_figure(p_acc, "Fig1_utility_boxplots_paired", width = 10, height = 5)
}

## ========== 3) FIGURE 2 — UMAP grids (paired layout + panel labels) ==========

png_paths <- list(
  ARF   = c(file.path(dir_arf, "umap_real_vs_arf.png"),
            file.path(dir_arf, "umap_arf.png")),
  CTGAN = c(file.path(dir_ctg, "umap_real_vs_ctgan.png"),
            file.path(dir_ctg, "umap_ctgan.png"))
)

# Build a single tile: load the PNG and add a panel label (a/b/c/d).
# No additional titles are drawn here (advisor requested removing extra titles).
build_umap_tile <- function(path, panel_letter) {
  if (!file.exists(path)) {
    return(
      ggdraw() +
        draw_label(glue("Missing: {basename(path)}"), size = 14)
    )
  }
  
  ggdraw() +
    draw_image(path, scale = 1) +
    # Panel label in upper-left corner
    draw_label(
      panel_letter,
      x = 0.02, y = 0.98,
      hjust = 0, vjust = 1,
      fontface = "bold",
      size = 18
    )
}

# Advisor request: place "Real vs CTGAN" as panel (b) in the upper-right.
# 2×2 layout:
# (a) upper-left  : Real vs ARF
# (b) upper-right : Real vs CTGAN
# (c) lower-left  : ARF-only
# (d) lower-right : CTGAN-only
p_umap <- plot_grid(
  build_umap_tile(png_paths$ARF[1],   "a"),
  build_umap_tile(png_paths$CTGAN[1], "b"),
  build_umap_tile(png_paths$ARF[2],   "c"),
  build_umap_tile(png_paths$CTGAN[2], "d"),
  ncol = 2,
  align = "hv"
)
save_figure(p_umap, "Fig2_umap_grids", width = 11, height = 8)
## ========== 4) FIGURE 3 — Privacy AUC (raw & normalized) ==========
auc_arf_raw <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf.csv"))
auc_arf_adj <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf_normalized.csv"))
auc_ctg_raw <- read_tbl(file.path(dir_ctg, "privacy_auc_cv_ctgan.csv"))
auc_ctg_adj <- read_tbl(file.path(dir_ctg, "privacy_auc_cv_ctgan_normalized.csv"))

if (!is.null(auc_arf_raw) & !is.null(auc_ctg_raw) & !is.null(auc_arf_adj) & !is.null(auc_ctg_adj)) {
  
  df_priv <- tibble(
    model = c("ARF","CTGAN"),
    AUC_raw = c(auc_arf_raw$.estimate[1], auc_ctg_raw$.estimate[1]),
    AUC_norm = c(auc_arf_adj$auc_report[1], auc_ctg_adj$auc_report[1])
  ) |>
    pivot_longer(cols = c(AUC_raw, AUC_norm), names_to = "type", values_to = "AUC") |>
    mutate(
      type = factor(type, levels = c("AUC_raw", "AUC_norm"))
    )
  
  dodge <- position_dodge(width = 0.6)
  
  p_priv <- ggplot(df_priv, aes(model, AUC, fill = type)) +
    geom_col(position = dodge, width = 0.55) +
    
    # Add AUC values on top of bars (advisor requested replacing the table)
    geom_text(
      aes(label = sprintf("%.3f", AUC)),
      position = dodge,
      vjust = -0.35,
      size = 4
    ) +
    
    geom_hline(yintercept = 0.5, linetype = 2) +
    
    # Add a bit of headroom so the "1.000" label is not clipped
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    
    scale_fill_manual(
      values = c("#9ecae1", "#6baed6"),
      labels = c("Raw AUC", "AUC normalized (≥0.5)")
    ) +
    
    labs(
      title = "Privacy: Domain-separator AUC (5×5 CV)",
      x = NULL, y = "AUC"
    ) +
    theme_pub +
    theme(legend.position = "top")
  
  save_figure(p_priv, "Fig3_privacy_auc", width = 7, height = 5)
}


## ========== 5) FIGURE 4 — Fidelity summary (BH-FDR %) ==========
fid_arf <- read_tbl(file.path(dir_arf, "fidelity_summary_arf.csv"))
fid_ctg <- read_tbl(file.path(dir_ctg, "fidelity_summary_ctgan.csv"))

if (!is.null(fid_arf) & !is.null(fid_ctg)) {
  df_fid <- bind_rows(
    fid_arf |> mutate(model="ARF"),
    fid_ctg |> mutate(model="CTGAN")
  ) |>
    select(model, ks_FDR_sig_pct, wil_FDR_sig_pct) |>
    pivot_longer(-model, names_to="test", values_to="pct") |>
    mutate(test = recode(test,
                         ks_FDR_sig_pct  = "KS (BH-FDR)",
                         wil_FDR_sig_pct = "Wilcoxon (BH-FDR)"))
  
  p_fid <- ggplot(df_fid, aes(test, pct, fill=model)) +
    geom_col(position=position_dodge(width=0.6), width=0.55) +
    scale_y_continuous(labels = percent_format(scale = 1), limits=c(0,100)) +
    labs(title="Fidelity: % genes significant (BH-FDR)", x=NULL, y="% significant") +
    theme_pub + theme(legend.position="top")
  save_figure(p_fid, "Fig4_fidelity_summary", width=8, height=5)
}

## ========== 6) FIGURE 5 — RF Importance Top-20 ==========
rf_real_arf <- read_tbl(file.path(dir_arf, "rf_importance_real.csv"))
rf_arf      <- read_tbl(file.path(dir_arf, "rf_importance_arf.csv"))
rf_real_ctg <- read_tbl(file.path(dir_ctg, "rf_importance_real.csv"))
rf_ctg      <- read_tbl(file.path(dir_ctg, "rf_importance_ctgan.csv"))

make_top20_bar <- function(tbl, title, value_col = "rf_importance") {
  if (is.null(tbl)) return(ggdraw() + draw_label(glue("Missing: {title}")))
  tbl %>%
    slice_head(n = min(20, nrow(.))) %>%
    mutate(gene = fct_reorder(gene, !!sym(value_col))) %>%
    ggplot(aes(gene, !!sym(value_col))) +
    geom_col() +
    coord_flip() +
    labs(title = title, x = "Gene", y = "Permutation Importance") +
    theme_pub +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}

# Clarify panel meaning: Real vs Synthetic within each generator.
# Remove the ambiguous word "run" from the Real panels.
p_arf_top <- plot_grid(
  make_top20_bar(rf_real_arf, "Real (RF importance)"),
  make_top20_bar(rf_arf,      "Synthetic (ARF)"),
  nrow = 1, rel_widths = c(1, 1)
)

p_ctg_top <- plot_grid(
  make_top20_bar(rf_real_ctg, "Real (RF importance)"),
  make_top20_bar(rf_ctg,      "Synthetic (CTGAN)"),
  nrow = 1, rel_widths = c(1, 1)
)

# Add a short header + optional note defining ARF/CTGAN in-figure (not required, but helpful)
p_rf <- plot_grid(
  ggdraw() + draw_label("Top-20 RF permutation importance", fontface = "bold", size = 16),
  p_arf_top,
  p_ctg_top,
  ncol = 1, rel_heights = c(0.12, 1, 1)
)

save_figure(p_rf, "Fig5_rf_importance_top20", width = 12, height = 10)

# Clarify panel meaning: Real vs Synthetic within each generator.
# Remove the ambiguous word "run" from the Real panels.
p_arf_top <- plot_grid(
  make_top20_bar(rf_real_arf, "Real (RF importance)"),
  make_top20_bar(rf_arf,      "Synthetic (ARF)"),
  nrow = 1, rel_widths = c(1, 1)
)

p_ctg_top <- plot_grid(
  make_top20_bar(rf_real_ctg, "Real (RF importance)"),
  make_top20_bar(rf_ctg,      "Synthetic (CTGAN)"),
  nrow = 1, rel_widths = c(1, 1)
)

# Add a short header + optional note defining ARF/CTGAN in-figure (not required, but helpful)
p_rf <- plot_grid(
  ggdraw() + draw_label("Top-20 RF permutation importance", fontface = "bold", size = 16),
  p_arf_top,
  p_ctg_top,
  ncol = 1, rel_heights = c(0.12, 1, 1)
)

save_figure(p_rf, "Fig5_rf_importance_top20", width = 12, height = 10)

## ========== 7) FIGURE 6 — NPDR Jaccard overlaps ==========
npdr_arf <- read_tbl(file.path(dir_arf, "npdr_pairwise_overlap_summary_real_arf_realplusarf.csv"))
npdr_ctg <- read_tbl(file.path(dir_ctg, "npdr_pairwise_overlap_summary_real_ctgan_realplusctgan.csv"))

plot_npdr_bar <- function(df, title){
  if (is.null(df)) return(ggdraw() + draw_label(glue("Missing: {title}")))
  df %>%
    mutate(pair = factor(pair, levels = pair)) %>%
    ggplot(aes(pair, jaccard)) +
    geom_col(width=0.6) +
    geom_text(aes(label=sprintf("n=%d", n_overlap)), vjust=-0.3, size=3.4) +
    scale_y_continuous(limits = c(0,1)) +
    labs(title = title, x=NULL, y="Jaccard") +
    theme_pub
}
p_npdr <- plot_grid(
  plot_npdr_bar(npdr_arf, "NPDR Top-K overlap — ARF"),
  plot_npdr_bar(npdr_ctg, "NPDR Top-K overlap — CTGAN"),
  nrow=1
)
save_figure(p_npdr, "Fig6_npdr_overlap_jaccard", width=10, height=4.5)

## ========== 8) FIGURE 7 — Dilution metric (rank drop on adding Synthetic) ==========
dil_arf  <- read_tbl(file.path(dir_arf, "dilution_detail_real_to_realplusARF_RF.csv")) |>
  mutate(model="ARF")
dil_ctg  <- read_tbl(file.path(dir_ctg, "dilution_detail_real_to_realplusCTGAN_RF.csv")) |>
  mutate(model="CTGAN")

if (!is.null(dil_arf) & !is.null(dil_ctg)) {
  dil_all <- bind_rows(dil_arf, dil_ctg)
  
  # Summary bars only (removed histogram panel per advisor feedback)
  p_dil2 <- dil_all |>
    group_by(model) |>
    summarise(mean_delta = mean(delta_rank, na.rm=TRUE),
              median_delta = median(delta_rank, na.rm=TRUE),
              pct_dropped  = mean(dropped==1, na.rm=TRUE)*100) |>
    pivot_longer(-model, names_to="metric", values_to="value") |>
    ggplot(aes(metric, value, fill=model)) +
    geom_col(position=position_dodge(width=0.7), width=0.6) +
    labs(title="Dilution summary", x=NULL, y="Value") +
    theme_pub +
    theme(axis.text.x = element_text(angle=20, hjust=1))
  
  save_figure(p_dil2, "Fig7_dilution_metric", width=10, height=5)
}
## ========== 9) FIGURE 8 — Stability (bootstrap Jaccard) ==========
stab_arf <- read_tbl(file.path(dir_arf, "jaccard_bootstrap_RF.csv")) |>
  mutate(model = "ARF")
stab_ctg <- read_tbl(file.path(dir_ctg, "jaccard_bootstrap_RF.csv")) |>
  mutate(model = "CTGAN")

if (!is.null(stab_arf) & !is.null(stab_ctg)) {
  stab_all <- bind_rows(stab_arf, stab_ctg) |>
    mutate(set = factor(set,
                        levels = c("Real_self","ARF_self","Real_vs_ARF","CTGAN_self","Real_vs_CTGAN")))
  
  p_stab <- ggplot(stab_all, aes(set, jaccard, fill = model)) +
    geom_violin(trim = FALSE, alpha = 0.9, position = position_dodge(width = 0.75)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.75)) +
    coord_cartesian(ylim = c(0, 0.3)) +  # <-- advisor request: reduce y-axis range
    scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
    labs(title = "Stability of Top-K (Jaccard over bootstraps)", x = NULL, y = "Jaccard") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  
  save_figure(p_stab, "Fig8_stability_jaccard", width = 10, height = 5.5)
}

message("All done. Figures saved to: ", fig_dir)

############################################################
# End of make-figures.R
############################################################



library(readr); library(dplyr)

arf_sum <- read_csv("~/Desktop/My paper/ARF_eval/table_acc_summary_arf.csv") %>%
  mutate(model = "ARF")
arf_t   <- read_csv("~/Desktop/My paper/ARF_eval/effectsize_real_vs_realplusarf.csv") %>%
  mutate(model = "ARF")

ctg_sum <- read_csv("~/Desktop/My paper/CTGAN_eval/table_acc_summary_ctgan.csv") %>%
  mutate(model = "CTGAN")
ctg_t   <- read_csv("~/Desktop/My paper/CTGAN_eval/effectsize_real_vs_realplusctgan.csv") %>%
  mutate(model = "CTGAN")

table1 <- bind_rows(
  arf_sum %>% left_join(arf_t, by="model"),
  ctg_sum %>% left_join(ctg_t,  by="model")
) %>%
  select(model, real_mean, real_sd, arf_mean = ctg_mean, arf_sd = ctg_sd,
         comb_mean, comb_sd, mean_diff, ci_low, ci_high, t, df, p_value, cohens_dz)

write_csv(table1, "~/Desktop/My paper/Table1_utility_summary_with_ttest.csv")
