############################################################
# regenerate_figures_colored.R
#
# Re-draws Fig1, Fidelity, Fig3, Fig5 (NPDR), Fig6, Fig7 with the
# shared ARF/CTGAN/TVAE color palette, top-right legends, and the
# wider "stretched" aspect ratio -- WITHOUT re-running any of the
# expensive RF/utility/stability computations.
#
# Reads the CSVs already saved by evaluate_all_generators.R under
# reviewer_revision/all_generators/{tables,details}/ and overwrites
# only the PNG files in .../figures/.
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})

## Shared theme: Times New Roman, horizontal axis labels, used in every figure.
FONT_FAMILY <- "Times New Roman"
THEME_BASE <- theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    strip.text = element_text(family = FONT_FAMILY, face = "bold", size = 12),
    legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12),
    legend.title = element_text(family = FONT_FAMILY, face = "bold", size = 12)
  )
## Always save with ragg so the Times New Roman font is embedded correctly.
save_fig <- function(filename, plot, width, height) {
  ggsave(filename, plot, width = width, height = height, dpi = 600, device = ragg::agg_png, bg = "white")
}

base_dir    <- path.expand("~/Desktop/My paper")
out_dir     <- file.path(base_dir, "reviewer_revision", "all_generators")
tables_dir  <- file.path(out_dir, "tables")
details_dir <- file.path(out_dir, "details")
figures_dir <- file.path(out_dir, "figures")

TOP_K <- 50

GENERATOR_COLORS <- c(
  ARF   = "#F8766D",
  CTGAN = "#00BFC4",
  TVAE  = "#E5A500"
)

## ---------- Fig1: Utility ----------
utility_wide <- read_csv(file.path(details_dir, "01_utility_replicates.csv"), show_col_types = FALSE)
utility_long <- utility_wide |>
  pivot_longer(cols = c(Real, Synthetic, Real_plus_Synthetic), names_to = "condition", values_to = "accuracy") |>
  mutate(
    condition = recode(condition, Real = "Real", Synthetic = "Synthetic", Real_plus_Synthetic = "Real+Synthetic"),
    condition = factor(condition, levels = c("Real", "Synthetic", "Real+Synthetic"))
  )
p_utility <- ggplot(utility_long, aes(x = condition, y = accuracy, fill = generator)) +
  geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 2.5, stroke = 1, color = "black") +
  facet_wrap(~ generator, nrow = 1) +
  scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
  THEME_BASE +
  labs(x = NULL, y = "RF OOB accuracy", title = NULL) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
        legend.position = "top", legend.justification = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.margin = margin(3, 8, 3, 8))
save_fig(file.path(figures_dir, "Fig1_Utility_all_generators.png"), p_utility, width = 13, height = 6.2)
cat("Saved Fig1_Utility_all_generators.png\n")

## ---------- Fig2: UMAP (single combined 2x3 panel figure) ----------
umap_path <- file.path(details_dir, "02_umap_coordinates.csv")
if (file.exists(umap_path)) {
  library(patchwork); library(cowplot)
  umap_all <- read_csv(umap_path, show_col_types = FALSE)
  ALL_UMAP_COLORS <- c(GENERATOR_COLORS, Real = "grey40")

  make_panel <- function(df, show_legend, extra_top = 0, extra_bottom = 0) {
    df$source <- factor(df$source, levels = c("Real", "ARF", "CTGAN", "TVAE"))
    df$class <- factor(df$class, levels = c("HC", "MDD"))
    ggplot(df, aes(UMAP1, UMAP2, color = source, shape = class)) +
      geom_point(alpha = 0.75, size = 1.8) +
      scale_color_manual(values = ALL_UMAP_COLORS, name = NULL, drop = FALSE, breaks = c("ARF", "CTGAN", "TVAE")) +
      scale_shape_manual(values = c(HC = 17, MDD = 15), name = NULL, drop = FALSE) +
      guides(
        color = guide_legend(order = 1, override.aes = list(shape = 15, size = 4)),
        shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3))
      ) +
      THEME_BASE +
      labs(x = "V1", y = "V2") +
      theme(panel.grid.minor = element_blank(),
            legend.box = "vertical",
            legend.position = if (show_legend) "top" else "none",
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
            legend.margin = margin(8, 14, 8, 14),
            legend.key.size = unit(0.65, "cm"),
            legend.spacing.x = unit(0.3, "cm"),
            legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12, margin = margin(r = 10)),
            plot.margin = margin(t = 5 + extra_top, b = 5 + extra_bottom, l = 5, r = 5))
  }

  top_row_data <- lapply(c("ARF", "CTGAN", "TVAE"), function(g) umap_all |> filter(panel == paste0("Real_vs_", g)))
  bottom_row_data <- lapply(c("ARF", "CTGAN", "TVAE"), function(g) umap_all |> filter(panel == paste0(g, "_only")) |> mutate(source = g))

  ## Extract ONE shared legend manually (avoids patchwork's guide-collection
  ## duplicating the legend across panels).
  legend_dummy <- expand.grid(
    UMAP1 = 0, UMAP2 = 0,
    source = c("ARF", "CTGAN", "TVAE"),
    class = c("HC", "MDD")
  )
  legend_source_plot <- make_panel(legend_dummy, show_legend = TRUE)
  shared_legend <- cowplot::get_legend(legend_source_plot)

  ## Arrange as 3 rows x 2 cols, grouped by generator, matching the original
  ## Fig.2 layout: (a) Real+ARF, (b) ARF only, (c) Real+CTGAN, (d) CTGAN only,
  ## (e) Real+TVAE, (f) TVAE only.
  panel_rows <- lapply(seq_along(c("ARF", "CTGAN", "TVAE")), function(i) {
    left  <- make_panel(top_row_data[[i]],    show_legend = FALSE, extra_bottom = 8)
    right <- make_panel(bottom_row_data[[i]], show_legend = FALSE, extra_bottom = 8)
    left | right
  })

  panels_grid <- (panel_rows[[1]] / panel_rows[[2]] / panel_rows[[3]]) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = 12),
          plot.tag.position = "top")

  legend_row <- cowplot::plot_grid(NULL, shared_legend, ncol = 2, rel_widths = c(0.55, 0.45))
  p_umap_combined <- cowplot::plot_grid(
    legend_row, panels_grid,
    ncol = 1, rel_heights = c(0.11, 1)
  )
  save_fig(file.path(figures_dir, "Fig2_UMAP_all_generators.png"), p_umap_combined, width = 12, height = 13.5)
  cat("Saved Fig2_UMAP_all_generators.png -- (a) Real+ARF (b) ARF only (c) Real+CTGAN (d) CTGAN only (e) Real+TVAE (f) TVAE only\n")
} else {
  cat("Skipped Fig2 (UMAP coordinates not found).\n")
}

## ---------- Fig4: RF Feature Importance ----------
reorder_within <- function(x, by, within, sep = "___") stats::reorder(paste(x, within, sep = sep), by)
scale_y_reordered <- function(...) scale_y_discrete(labels = function(x) gsub("___.*$", "", x), ...)
TOP_BAR <- 20
rf_real_path <- file.path(details_dir, "05_rf_importance_Real.csv")
if (file.exists(rf_real_path)) {
  rf_datasets <- c("Real", "ARF", "CTGAN", "TVAE")
  RF_PLOT_COLORS <- c(GENERATOR_COLORS, Real = "grey40")
  rf_data <- lapply(rf_datasets, function(nm) {
    read_csv(file.path(details_dir, paste0("05_rf_importance_", nm, ".csv")), show_col_types = FALSE) |>
      slice_head(n = TOP_BAR) |> mutate(dataset = nm)
  })
  names(rf_data) <- rf_datasets

  make_rf_panel <- function(df, show_legend = FALSE) {
    df <- df |> mutate(gene_facet = reorder(gene, importance),
                        dataset = factor(dataset, levels = rf_datasets))
    ggplot(df, aes(importance, gene_facet, fill = dataset)) +
      geom_col() +
      scale_fill_manual(values = RF_PLOT_COLORS, name = NULL, drop = FALSE) +
      THEME_BASE + labs(x = "Permutation importance", y = "Gene") +
      theme(panel.grid.minor = element_blank(),
            legend.position = if (show_legend) "top" else "none",
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
            legend.margin = margin(8, 18, 8, 18),
            legend.key.size = unit(0.65, "cm"),
            legend.spacing.x = unit(0.4, "cm"),
            legend.text = element_text(family = FONT_FAMILY, face = "bold", size = 12, margin = margin(r = 12)))
  }

  ## Shared legend built from a dummy dataset containing all 4 categories,
  ## so every color swatch renders (a single real panel may be missing a
  ## category's data and silently drop its swatch -- see the Fig.2 fix).
  legend_dummy <- tibble(gene = paste0("g", seq_along(rf_datasets)),
                          importance = 1, dataset = rf_datasets)
  shared_legend <- cowplot::get_legend(make_rf_panel(legend_dummy, show_legend = TRUE))

  panels <- lapply(rf_datasets, function(nm) make_rf_panel(rf_data[[nm]]))
  rf_grid <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = 12),
          plot.tag.position = "top")
  legend_row <- cowplot::plot_grid(NULL, shared_legend, ncol = 2, rel_widths = c(0.35, 0.65))
  p_rf <- cowplot::plot_grid(legend_row, rf_grid, ncol = 1, rel_heights = c(0.08, 1))
  save_fig(file.path(figures_dir, "Fig4_RF_importance_top20_all_generators.png"), p_rf, width = 13, height = 12)
  cat("Saved Fig4_RF_importance_top20_all_generators.png -- (a) Real (b) ARF (c) CTGAN (d) TVAE\n")
} else {
  cat("Skipped Fig4 (RF importance detail files not found).\n")
}

## ---------- Fig3: Privacy (Normalized domain-separator AUC only) ----------
privacy_summary <- read_csv(file.path(tables_dir, "04_privacy_summary.csv"), show_col_types = FALSE)
privacy_plot_df <- privacy_summary |>
  mutate(generator = factor(generator, levels = c("ARF", "CTGAN", "TVAE")))
p_privacy <- ggplot(privacy_plot_df, aes(generator, normalized_auc, fill = generator)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.3f", normalized_auc)),
            vjust = -0.5, size = 4, family = FONT_FAMILY, fontface = "bold") +
  scale_fill_manual(values = GENERATOR_COLORS, guide = "none") +
  coord_cartesian(ylim = c(0, 1.08)) + THEME_BASE +
  labs(x = NULL, y = "Normalized domain-separator AUC", title = NULL) +
  theme(panel.grid.minor = element_blank())
save_fig(file.path(figures_dir, "Fig3_Privacy_all_generators.png"), p_privacy, width = 6.5, height = 5.5)
cat("Saved Fig3_Privacy_all_generators.png\n")

## ---------- Fig5: NPDR (only if the file exists -- may be missing if npdr wasn't run) ----------
npdr_path <- file.path(tables_dir, "06_npdr_feature_alignment_summary.csv")
if (file.exists(npdr_path)) {
  npdr_alignment_summary <- read_csv(npdr_path, show_col_types = FALSE) |>
    mutate(comparison = factor(comparison, levels = c(
      "Real vs Synthetic", "Real vs Real+Synthetic", "Synthetic vs Real+Synthetic"
    )))
  p_npdr <- ggplot(npdr_alignment_summary, aes(comparison, jaccard, fill = generator)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) + THEME_BASE +
    scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
    labs(x = NULL, y = paste0("NPDR Top-", TOP_K, " Jaccard"), fill = "Generator", title = NULL) +
    theme(panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "right")
  save_fig(file.path(figures_dir, "Fig5_NPDR_alignment_all_generators.png"), p_npdr, width = 9, height = 5)
  cat("Saved Fig5_NPDR_alignment_all_generators.png\n")
} else {
  cat("Skipped Fig5 (NPDR summary not found).\n")
}

## ---------- Fig6: Dilution ----------
ranking_dilution_summary <- read_csv(file.path(tables_dir, "07_ranking_dilution_summary_K50.csv"), show_col_types = FALSE)
dilution_plot_df <- ranking_dilution_summary |>
  transmute(generator = gsub("^Real vs Real\\+", "", Comparison),
            `Mean absolute rank shift` = Mean_absolute_rank_shift,
            `Median absolute rank shift` = Median_absolute_rank_shift,
            `Top-50 dropout (%)` = Top50_dropout_pct) |>
  pivot_longer(-generator, names_to = "metric", values_to = "value") |>
  mutate(generator = factor(generator, levels = c("ARF", "CTGAN", "TVAE")),
         metric = factor(metric, levels = c("Mean absolute rank shift", "Median absolute rank shift", "Top-50 dropout (%)")))
p_dilution <- ggplot(dilution_plot_df, aes(generator, value, fill = generator)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.5, size = 4, family = FONT_FAMILY, fontface = "bold") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = GENERATOR_COLORS, name = NULL) +
  THEME_BASE +
  labs(x = NULL, y = "Value", title = NULL) +
  theme(panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top", legend.justification = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.margin = margin(3, 8, 3, 8))
save_fig(file.path(figures_dir, "Fig6_Dilution_all_generators.png"), p_dilution, width = 12, height = 5)
cat("Saved Fig6_Dilution_all_generators.png\n")

## ---------- Fig7: Stability (violin, matching original style) ----------
gen_names <- c("ARF", "CTGAN", "TVAE")
stability_replicates <- read_csv(file.path(details_dir, "08_stability_replicates.csv"), show_col_types = FALSE) |>
  mutate(fill_group = case_when(
    comparison == "Real" ~ "Real",
    comparison %in% gen_names ~ comparison,
    TRUE ~ sub("^Real vs ", "", comparison)
  ))
STABILITY_COLORS <- c(GENERATOR_COLORS, Real = "grey50")
STABILITY_ORDER <- c("Real", "ARF", "Real vs ARF", "CTGAN", "Real vs CTGAN", "TVAE", "Real vs TVAE")
stability_replicates <- stability_replicates |>
  mutate(comparison = factor(comparison, levels = STABILITY_ORDER),
         fill_group = factor(fill_group, levels = c("Real", gen_names)))
p_stability <- ggplot(stability_replicates, aes(comparison, jaccard, fill = fill_group)) +
  geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = STABILITY_COLORS, name = NULL, breaks = gen_names) +
  THEME_BASE +
  labs(x = NULL, y = paste0("Bootstrap Top-", TOP_K, " Jaccard"), title = NULL) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        legend.position = "top", legend.justification = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.margin = margin(3, 8, 3, 8))
save_fig(file.path(figures_dir, "Fig7_Stability_all_generators.png"), p_stability, width = 12, height = 6)
cat("Saved Fig7_Stability_all_generators.png\n")

cat("\nAll figures regenerated with the shared color palette (no RF re-computation needed).\n")
