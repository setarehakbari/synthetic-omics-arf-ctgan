############################################################
# make-figures.R
# Output: ~/Desktop/My paper/results
############################################################

## ========== 0) Setup ==========
base_dir <- path.expand("~/Desktop/My paper")
dir_arf  <- file.path(base_dir, "ARF_eval")
dir_ctg  <- file.path(base_dir, "CTGAN_eval")
res_dir  <- file.path(base_dir, "results")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

## ========== 1) Packages ==========
need <- c("tidyverse", "cowplot", "scales", "glue", "patchwork", "magick")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Please install package: '%s'", p))
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(scales)
  library(glue)
  library(patchwork)
  library(magick)
  library(grid)
})

has_gt      <- requireNamespace("gt", quietly = TRUE)
has_webshot <- requireNamespace("webshot2", quietly = TRUE)

base_family <- "Times"

## ========== 2) Helpers ==========
read_tbl <- function(path) {
  if (file.exists(path)) {
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  } else {
    message("Missing: ", path)
    NULL
  }
}

save_figure <- function(plot_obj, name, width = 8, height = 5, dpi = 300) {
  png_path <- file.path(res_dir, paste0(name, ".png"))
  pdf_path <- file.path(res_dir, paste0(name, ".pdf"))
  
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi)
  
  ok <- TRUE
  tryCatch(
    {
      ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf)
    },
    error = function(e) {
      ok <<- FALSE
    }
  )
  if (!ok) {
    ggsave(pdf_path, plot_obj, width = width, height = height, device = "pdf")
  }
}

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

############################################################
## FIGURE 1 — Utility (RF OOB accuracy)
############################################################
acc_arf <- read_tbl(file.path(dir_arf, "rf_accuracies_arf_vs_real_combined.csv")) |>
  mutate(model = "ARF") |>
  rename(synth = arf, combined = combined_arf)

acc_ctg <- read_tbl(file.path(dir_ctg, "rf_accuracies_ctgan_vs_real_combined.csv")) |>
  mutate(model = "CTGAN") |>
  rename(synth = ctgan, combined = combined_ctgan)

if (!is.null(acc_arf) && !is.null(acc_ctg)) {
  
  acc_long <- bind_rows(acc_arf, acc_ctg) |>
    pivot_longer(c(real, synth, combined), names_to = "set", values_to = "accuracy") |>
    mutate(
      set = factor(
        set,
        levels = c("real", "synth", "combined"),
        labels = c("Real", "Synthetic", "Real+Synthetic")
      )
    )
  
  p_acc_base <- ggplot(acc_long, aes(set, accuracy, fill = model)) +
    geom_violin(trim = FALSE, alpha = 0.9) +
    geom_boxplot(width = 0.16, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3, stroke = 1.1) +
    facet_wrap(~model, nrow = 1) +
    coord_cartesian(ylim = c(0.45, 0.70)) +
    scale_fill_manual(values = c("ARF" = "#F8766D", "CTGAN" = "#00BFC4")) +
    labs(x = NULL, y = "Accuracy (OOB)") +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_text(family = base_family, size = 15, face = "plain"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  p_acc <- ggdraw() +
    draw_plot(p_acc_base, x = 0, y = 0, width = 1, height = 0.92) +
    draw_grob(
      rectGrob(gp = gpar(fill = "white", col = "black", lwd = 0.8)),
      x = 0.76, y = 0.93, width = 0.23, height = 0.065
    ) +
    draw_grob(
      rectGrob(gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)),
      x = 0.80, y = 0.947, width = 0.014, height = 0.020
    ) +
    draw_text(
      "ARF",
      x = 0.845, y = 0.957,
      family = base_family, size = 15, fontface = "plain", color = "black"
    ) +
    draw_grob(
      rectGrob(gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)),
      x = 0.885, y = 0.947, width = 0.014, height = 0.020
    ) +
    draw_text(
      "CTGAN",
      x = 0.950, y = 0.957,
      family = base_family, size = 15, fontface = "plain", color = "black"
    )
  
  save_figure(p_acc, "Fig1_utility_boxplots", 10, 5.4)
}
############################################################
## FIGURE 2 — UMAP grids from existing PNG panels
############################################################

library(magick)
library(cowplot)
library(ggplot2)
library(glue)
library(grid)

base_family <- "Times"

png_paths <- list(
  ARF = c(
    file.path(dir_arf, "umap_real_vs_arf.png"),
    file.path(dir_arf, "umap_arf.png")
  ),
  CTGAN = c(
    file.path(dir_ctg, "umap_real_vs_ctgan.png"),
    file.path(dir_ctg, "umap_ctgan.png")
  )
)

# ----------------------------------------------------------
# Crop helper
# ----------------------------------------------------------
crop_umap_image <- function(path,
                            crop_left = 42,
                            crop_right = 170,
                            crop_top = 72,
                            crop_bottom = 42) {
  img  <- magick::image_read(path)
  info <- magick::image_info(img)
  
  new_w <- info$width  - crop_left - crop_right
  new_h <- info$height - crop_top  - crop_bottom
  
  magick::image_crop(
    img,
    magick::geometry_area(
      width  = new_w,
      height = new_h,
      x_off  = crop_left,
      y_off  = crop_top
    )
  )
}

# ----------------------------------------------------------
# Build one tile
# ----------------------------------------------------------
build_umap_tile <- function(path, panel_label,
                            crop_left = 42,
                            crop_right = 170,
                            crop_top = 72,
                            crop_bottom = 42,
                            img_x = 0.13,
                            img_y = 0.13,
                            img_width = 0.80,
                            img_height = 0.78,
                            add_axes = TRUE) {
  if (!file.exists(path)) {
    return(
      ggdraw() +
        draw_label(
          glue("Missing: {basename(path)}"),
          x = 0.5, y = 0.5,
          fontfamily = base_family,
          size = 15,
          fontface = "plain"
        )
    )
  }
  
  img <- crop_umap_image(
    path,
    crop_left   = crop_left,
    crop_right  = crop_right,
    crop_top    = crop_top,
    crop_bottom = crop_bottom
  )
  
  p <- ggdraw() +
    draw_label(
      panel_label,
      x = 0.5, y = 0.985,
      hjust = 0.5, vjust = 1,
      fontfamily = base_family,
      size = 15,
      fontface = "plain"
    ) +
    draw_image(
      img,
      x = img_x, y = img_y,
      width = img_width, height = img_height,
      scale = 1
    )
  
  if (add_axes) {
    p <- p +
      draw_label(
        "V1",
        x = 0.53, y = 0.055,
        hjust = 0.5, vjust = 0.5,
        fontfamily = base_family,
        size = 15,
        fontface = "plain"
      ) +
      draw_label(
        "V2",
        x = 0.045, y = 0.50,
        angle = 90,
        hjust = 0.5, vjust = 0.5,
        fontfamily = base_family,
        size = 15,
        fontface = "plain"
      )
  }
  
  p
}

# ----------------------------------------------------------
# Row 1
# ----------------------------------------------------------
row_1 <- plot_grid(
  build_umap_tile(
    png_paths$ARF[1],
    "(a)",
    crop_left = 42,
    crop_right = 170,
    crop_top = 72,
    crop_bottom = 42,
    img_x = 0.13,
    img_y = 0.13,
    img_width = 0.80,
    img_height = 0.78,
    add_axes = TRUE
  ),
  build_umap_tile(
    png_paths$ARF[2],
    "(b)",
    crop_left = 36,
    crop_right = 170,
    crop_top = 72,
    crop_bottom = 28,
    img_x = 0.13,
    img_y = 0.15,
    img_width = 0.80,
    img_height = 0.78,
    add_axes = TRUE
  ),
  nrow = 1,
  rel_widths = c(1, 1)
)

# ----------------------------------------------------------
# Row 2
# ----------------------------------------------------------
row_2 <- plot_grid(
  build_umap_tile(
    png_paths$CTGAN[1],
    "(c)",
    crop_left = 42,
    crop_right = 260,
    crop_top = 78,     
    crop_bottom = 58,   
    img_x = 0.13,
    img_y = 0.14,        
    img_width = 0.80,
    img_height = 0.76,  
    add_axes = TRUE
  ),
  build_umap_tile(
    png_paths$CTGAN[2],
    "(d)",
    crop_left = 42,
    crop_right = 170,
    crop_top = 72,
    crop_bottom = 42,
    img_x = 0.13,
    img_y = 0.13,
    img_width = 0.80,
    img_height = 0.78,
    add_axes = TRUE
  ),
  nrow = 1,
  rel_widths = c(1, 1)
)

main_umap <- plot_grid(
  row_1,
  row_2,
  ncol = 1,
  rel_heights = c(1, 1)
)

# ----------------------------------------------------------
# Final figure with shared legend
# ----------------------------------------------------------
p_umap <- ggdraw() +
  draw_plot(main_umap, x = 0, y = 0, width = 1, height = 0.93) +
  
  draw_grob(
    rectGrob(
      gp = gpar(fill = "white", col = "black", lwd = 0.8)
    ),
    x = 0.70, y = 0.935, width = 0.28, height = 0.065
  ) +
  
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)
    ),
    x = 0.72, y = 0.959, width = 0.012, height = 0.016
  ) +
  draw_label(
    "ARF",
    x = 0.755, y = 0.969,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)
    ),
    x = 0.80, y = 0.959, width = 0.012, height = 0.016
  ) +
  draw_label(
    "CTGAN",
    x = 0.865, y = 0.969,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  
  draw_label(
    "▲",
    x = 0.73, y = 0.948,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  draw_label(
    "HC",
    x = 0.760, y = 0.948,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  
  draw_label(
    "■",
    x = 0.81, y = 0.948,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  draw_label(
    "MDD",
    x = 0.850, y = 0.948,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  )

save_figure(p_umap, "Fig2_umap_grids", 12.5, 9.3)

############################################################
## FIGURE 3 — Privacy AUC
############################################################
auc_arf_raw <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf.csv"))
auc_arf_adj <- read_tbl(file.path(dir_arf, "privacy_auc_cv_arf_normalized.csv"))
auc_ctg_raw <- read_tbl(file.path(dir_ctg, "privacy_auc_cv_ctgan.csv"))
auc_ctg_adj <- read_tbl(file.path(dir_ctg, "privacy_auc_cv_ctgan_normalized.csv"))

if (!is.null(auc_arf_raw) && !is.null(auc_ctg_raw) &&
    !is.null(auc_arf_adj) && !is.null(auc_ctg_adj)) {
  
  df_priv <- tibble(
    model = c("ARF", "CTGAN"),
    AUC_raw  = c(auc_arf_raw$.estimate[1], auc_ctg_raw$.estimate[1]),
    AUC_norm = c(auc_arf_adj$auc_report[1], auc_ctg_adj$auc_report[1])
  ) |>
    pivot_longer(c(AUC_raw, AUC_norm), names_to = "type", values_to = "AUC") |>
    mutate(
      type = factor(
        type,
        levels = c("AUC_raw", "AUC_norm"),
        labels = c("Raw AUC", "Normalized AUC")
      ),
      model = factor(model, levels = c("ARF", "CTGAN")),
      label = sprintf("%.2f", AUC)
    )
  
  dodge_width <- 0.78
  
  p_priv_base <- ggplot(df_priv, aes(type, AUC, fill = model)) +
    geom_col(
      position = position_dodge(width = dodge_width),
      width = 0.34,
      color = "black",
      linewidth = 0.4
    ) +
    geom_text(
      aes(label = label),
      position = position_dodge(width = dodge_width),
      vjust = -0.35,
      family = base_family,
      size = 4.5
    ) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(
      limits = c(0, 1.08),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(values = c("ARF" = "#F8766D", "CTGAN" = "#00BFC4")) +
    labs(x = NULL, y = "AUC") +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_text(family = base_family, size = 15, face = "plain"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  p_priv <- ggdraw() +
    draw_plot(p_priv_base, x = 0, y = 0, width = 1, height = 0.92) +
    draw_grob(
      rectGrob(gp = gpar(fill = "white", col = "black", lwd = 0.8)),
      x = 0.72, y = 0.93, width = 0.23, height = 0.065
    ) +
    draw_grob(
      rectGrob(gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)),
      x = 0.74, y = 0.952, width = 0.014, height = 0.018
    ) +
    draw_label(
      "ARF",
      x = 0.790, y = 0.962,
      fontfamily = base_family, size = 15, fontface = "plain"
    ) +
    draw_grob(
      rectGrob(gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)),
      x = 0.83, y = 0.952, width = 0.014, height = 0.018
    ) +
    draw_label(
      "CTGAN",
      x = 0.900, y = 0.962,
      fontfamily = base_family, size = 15, fontface = "plain"
    )
  
  save_figure(p_priv, "Fig3_privacy_auc", 7.6, 5.4)
}

## =========================================================
## FIGURE 4 — RF Top-20
## =========================================================

library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(cowplot)
library(glue)
library(rlang)
library(grid)

# For macOS
base_family <- "Times"

rf_real_arf <- read_tbl(file.path(dir_arf, "rf_importance_real.csv"))
rf_arf      <- read_tbl(file.path(dir_arf, "rf_importance_arf.csv"))
rf_real_ctg <- read_tbl(file.path(dir_ctg, "rf_importance_real.csv"))
rf_ctg      <- read_tbl(file.path(dir_ctg, "rf_importance_ctgan.csv"))

make_top20 <- function(tbl, fill_color, value_col = "rf_importance") {
  if (is.null(tbl)) {
    return(
      ggdraw() +
        draw_label(
          "Missing panel",
          fontfamily = base_family,
          fontface = "plain",
          size = 15
        )
    )
  }
  
  tbl |>
    slice_head(n = min(20, nrow(tbl))) |>
    mutate(gene = fct_reorder(gene, !!sym(value_col))) |>
    ggplot(aes(gene, !!sym(value_col))) +
    geom_col(fill = fill_color) +
    coord_flip() +
    labs(
      x = "Gene",
      y = "Permutation Importance"
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_blank(),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      legend.title = element_text(family = base_family, size = 15, face = "plain"),
      legend.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

add_panel_label <- function(p, label) {
  ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_label(
      label,
      x = 0.5, y = 0.99,
      hjust = 0.5, vjust = 1,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    )
}

# ARF panels = red
p_a <- add_panel_label(
  make_top20(rf_real_arf, fill_color = "#F8766D"),
  "(a)"
)

p_b <- add_panel_label(
  make_top20(rf_arf, fill_color = "#F8766D"),
  "(b)"
)

# CTGAN panels = blue
p_c <- add_panel_label(
  make_top20(rf_real_ctg, fill_color = "#00BFC4"),
  "(c)"
)

p_d <- add_panel_label(
  make_top20(rf_ctg, fill_color = "#00BFC4"),
  "(d)"
)

# Base 2x2 panel layout
p_rf_base <- plot_grid(
  p_a, p_b,
  p_c, p_d,
  ncol = 2,
  align = "hv"
)

# Final figure with global legend in the top-right
p_rf <- ggdraw() +
  # lower the main plot slightly to create a top band for legend
  draw_plot(p_rf_base, x = 0, y = 0, width = 1, height = 0.92) +
  
  # legend box (larger, top-right)
  draw_grob(
    rectGrob(
      gp = gpar(fill = "white", col = "black", lwd = 0.8)
    ),
    x = 0.73, y = 0.93, width = 0.24, height = 0.07
  ) +
  
  # ARF color square
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)
    ),
    x = 0.75, y = 0.952, width = 0.014, height = 0.018
  ) +
  draw_label(
    "ARF",
    x = 0.800, y = 0.962,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  
  # CTGAN color square
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)
    ),
    x = 0.84, y = 0.952, width = 0.014, height = 0.018
  ) +
  draw_label(
    "CTGAN",
    x = 0.910, y = 0.962,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  )

save_figure(p_rf, "Fig4_rf_importance_top20", 12, 10.4)

## =========================================================
## FIGURE 5 — NPDR overlaps
## =========================================================

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(glue)
library(grid)
library(stringr)

# For macOS
base_family <- "Times"

npdr_arf <- read_tbl(file.path(dir_arf, "npdr_pairwise_overlap_summary_real_arf_realplusarf.csv"))
npdr_ctg <- read_tbl(file.path(dir_ctg, "npdr_pairwise_overlap_summary_real_ctgan_realplusctgan.csv"))

plot_npdr_bar <- function(df, fill_color) {
  if (is.null(df)) {
    return(
      ggdraw() +
        draw_label(
          "Missing panel",
          fontfamily = base_family,
          fontface = "plain",
          size = 15
        )
    )
  }
  
  df_plot <- df |>
    mutate(
      pair = as.character(pair),
      pair = gsub("_", " ", pair),
      pair = gsub("\\bvs\\b", "against", pair, ignore.case = TRUE),
      pair = str_wrap(pair, width = 16),
      pair = factor(pair, levels = unique(pair))
    )
  
  ggplot(df_plot, aes(pair, jaccard)) +
    geom_col(width = 0.28, fill = fill_color) +
    geom_text(
      aes(label = sprintf("%.2f", jaccard)),
      vjust = -0.45,
      family = base_family,
      size = 4.8
    ) +
    scale_x_discrete(
      expand = expansion(add = c(1.2, 1.2))
    ) +
    scale_y_continuous(
      limits = c(0, 1.10),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Jaccard"
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      legend.title = element_text(family = base_family, size = 15, face = "plain"),
      legend.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_blank(),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        family = base_family,
        size = 15,
        face = "plain"
      ),
      plot.margin = margin(t = 18, r = 12, b = 12, l = 12)
    )
}

add_panel_label <- function(p, label) {
  ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_label(
      label,
      x = 0.5, y = 0.99,
      hjust = 0.5, vjust = 1,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    )
}

# (a) ARF = red (top)
p_a <- add_panel_label(
  plot_npdr_bar(npdr_arf, fill_color = "#F8766D"),
  "(a)"
)

# (b) CTGAN = blue (bottom)
p_b <- add_panel_label(
  plot_npdr_bar(npdr_ctg, fill_color = "#00BFC4"),
  "(b)"
)

# stack panels vertically
p_npdr_base <- plot_grid(
  p_a,
  p_b,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v"
)

# final figure with shared legend at the top
p_npdr <- ggdraw() +
  draw_plot(p_npdr_base, x = 0, y = 0, width = 1, height = 0.88) +
  
  # legend box
  draw_grob(
    rectGrob(
      gp = gpar(fill = "white", col = "black", lwd = 0.8)
    ),
    x = 0.67, y = 0.920, width = 0.28, height = 0.075
  ) +
  
  # ARF square
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)
    ),
    x = 0.70, y = 0.943, width = 0.014, height = 0.018
  ) +
  draw_label(
    "ARF",
    x = 0.755, y = 0.953,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  ) +
  
  # CTGAN square
  draw_grob(
    rectGrob(
      gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)
    ),
    x = 0.80, y = 0.943, width = 0.014, height = 0.018
  ) +
  draw_label(
    "CTGAN",
    x = 0.875, y = 0.953,
    fontfamily = base_family,
    fontface = "plain",
    size = 15
  )

save_figure(p_npdr, "Fig5_npdr_overlap_jaccard", 9.2, 9.6)

## =========================================================
## FIGURE 6 — Dilution summary
## =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)
library(grid)

# For macOS
base_family <- "Times"

dil_arf <- read_tbl(file.path(dir_arf, "dilution_detail_real_to_realplusARF_RF.csv")) |>
  mutate(model = "ARF")

dil_ctg <- read_tbl(file.path(dir_ctg, "dilution_detail_real_to_realplusCTGAN_RF.csv")) |>
  mutate(model = "CTGAN")

if (!is.null(dil_arf) && !is.null(dil_ctg)) {
  
  dil_all <- bind_rows(dil_arf, dil_ctg)
  
  p_d_data <- dil_all |>
    group_by(model) |>
    summarise(
      mean_delta   = mean(delta_rank, na.rm = TRUE),
      median_delta = median(delta_rank, na.rm = TRUE),
      pct_dropped  = mean(dropped == 1, na.rm = TRUE) * 100,
      .groups = "drop"
    ) |>
    pivot_longer(-model, names_to = "metric", values_to = "value") |>
    mutate(
      metric = factor(
        metric,
        levels = c("mean_delta", "median_delta", "pct_dropped"),
        labels = c("Mean rank shift", "Median rank shift", "Top-K dropout (%)")
      )
    )
  
  p_d_base <- ggplot(p_d_data, aes(metric, value, fill = model)) +
    geom_col(position = position_dodge(0.72), width = 0.60) +
    geom_text(
      aes(label = sub("\\.0$", "", sprintf("%.1f", value))),
      position = position_dodge(0.72),
      vjust = -0.35,
      family = base_family,
      size = 5
    ) +
    scale_fill_manual(values = c("ARF" = "#F8766D", "CTGAN" = "#00BFC4")) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.14))
    ) +
    labs(
      x = NULL,
      y = "Value"
    ) +
    coord_cartesian(clip = "off") +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      legend.title = element_text(family = base_family, size = 15, face = "plain"),
      legend.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_blank(),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        family = base_family,
        size = 15,
        face = "plain"
      ),
      plot.margin = margin(t = 22, r = 12, b = 12, l = 12)
    )
  
  p_d <- ggdraw() +
    draw_plot(p_d_base, x = 0, y = 0, width = 1, height = 0.90) +
    
    # legend box
    draw_grob(
      rectGrob(
        gp = gpar(fill = "white", col = "black", lwd = 0.8)
      ),
      x = 0.72, y = 0.92, width = 0.24, height = 0.075
    ) +
    
    # ARF square
    draw_grob(
      rectGrob(
        gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)
      ),
      x = 0.74, y = 0.945, width = 0.014, height = 0.018
    ) +
    draw_label(
      "ARF",
      x = 0.79, y = 0.955,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    ) +
    
    # CTGAN square
    draw_grob(
      rectGrob(
        gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)
      ),
      x = 0.83, y = 0.945, width = 0.014, height = 0.018
    ) +
    draw_label(
      "CTGAN",
      x = 0.90, y = 0.955,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    )
  
  save_figure(p_d, "Fig6_dilution_metric", 10.2, 5.8)
}
## =========================================================
## FIGURE 7 — Stability (bootstrap Jaccard)
## =========================================================

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(grid)

# For macOS
base_family <- "Times"

stab_arf <- read_tbl(file.path(dir_arf, "jaccard_bootstrap_RF.csv")) |>
  mutate(model = "ARF")

stab_ctg <- read_tbl(file.path(dir_ctg, "jaccard_bootstrap_RF.csv")) |>
  mutate(model = "CTGAN")

if (!is.null(stab_arf) && !is.null(stab_ctg)) {
  
  stab_all <- bind_rows(stab_arf, stab_ctg) |>
    mutate(
      set = factor(
        set,
        levels = c("Real_self", "ARF_self", "Real_vs_ARF", "CTGAN_self", "Real_vs_CTGAN"),
        labels = c("Real", "ARF", "Real against ARF", "CTGAN", "Real against CTGAN")
      )
    )
  
  p_stab_base <- ggplot(stab_all, aes(set, jaccard, fill = model)) +
    geom_violin(
      trim = FALSE,
      alpha = 0.9,
      position = position_dodge(width = 0.8)
    ) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.35, by = 0.05)
    ) +
    coord_cartesian(ylim = c(0, 0.35))+
    scale_fill_manual(values = c("ARF" = "#F8766D", "CTGAN" = "#00BFC4")) +
    labs(
      x = NULL,
      y = "Jaccard"
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      text = element_text(family = base_family, size = 15, face = "plain"),
      axis.title = element_text(family = base_family, size = 15, face = "plain"),
      axis.text = element_text(family = base_family, size = 15, face = "plain"),
      legend.title = element_text(family = base_family, size = 15, face = "plain"),
      legend.text = element_text(family = base_family, size = 15, face = "plain"),
      plot.title = element_blank(),
      strip.text = element_text(family = base_family, size = 15, face = "plain"),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        family = base_family,
        size = 15,
        face = "plain"
      ),
      plot.margin = margin(t = 18, r = 12, b = 12, l = 12)
    )
  
  p_stab <- ggdraw() +
    draw_plot(p_stab_base, x = 0, y = 0, width = 1, height = 0.90) +
    
    # legend box
    draw_grob(
      rectGrob(
        gp = gpar(fill = "white", col = "black", lwd = 0.8)
      ),
      x = 0.70, y = 0.915, width = 0.26, height = 0.075
    ) +
    
    # ARF square
    draw_grob(
      rectGrob(
        gp = gpar(fill = "#F8766D", col = "black", lwd = 0.6)
      ),
      x = 0.72, y = 0.938, width = 0.014, height = 0.018
    ) +
    draw_label(
      "ARF",
      x = 0.775, y = 0.948,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    ) +
    
    # CTGAN square
    draw_grob(
      rectGrob(
        gp = gpar(fill = "#00BFC4", col = "black", lwd = 0.6)
      ),
      x = 0.82, y = 0.938, width = 0.014, height = 0.018
    ) +
    draw_label(
      "CTGAN",
      x = 0.895, y = 0.948,
      fontfamily = base_family,
      fontface = "plain",
      size = 15
    )
  
  save_figure(p_stab, "Fig7_stability_jaccard", 10.6, 5.8)
}

#####################################################
message("All done. Results saved to: ", res_dir)




