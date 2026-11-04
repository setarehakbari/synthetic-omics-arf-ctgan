## ================================================================
## Fig.1 -- Overview of the proposed evaluation framework (workflow diagram)
## Reviewer 2, Major Comment 2
## ================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

## ---- EDIT THIS to your output folder ----
base_dir    <- Sys.getenv("SYNTH_EVAL_BASE_DIR", unset = ".")
figures_dir <- file.path(base_dir, "reviewer_revision", "all_generators", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

FONT_FAMILY <- "Times New Roman"
FONT_PT     <- 16   ## intentionally larger than the 12pt used in the data figures --
                     ## this is a schematic diagram meant to be readable at a glance
TXT_SIZE    <- FONT_PT / .pt   ## ggplot2 geom_text size is in mm; .pt converts pt -> mm units

GENERATOR_COLORS <- c(ARF = "#F8766D", CTGAN = "#00BFC4", TVAE = "#E5A500")
BOOKEND_FILL   <- "#D6EAF8"   ## real-data / comparative-analysis boxes
BOOKEND_BORDER <- "#1B4F72"
FRAMEWORK_FILL <- "#EAF2F8"
FRAMEWORK_BORDER <- "#1B4F72"
METRIC_FILL    <- "#AED6F1"
METRIC_BORDER  <- "#1B4F72"
GREY_TXT       <- "#1B4F72"

## ---------- helper: rectangle + centered label ----------
rects  <- list()
labels <- list()
add_box <- function(x, y, w, h, text, fill = "white", color = "black", lw = 1.1,
                     textcolor = "black") {
  rects[[length(rects) + 1]] <<- data.frame(
    xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
    fill = fill, color = color, lw = lw
  )
  labels[[length(labels) + 1]] <<- data.frame(x = x, y = y, label = text, textcolor = textcolor)
}

## ---------- Stage 1: Real data ----------
add_box(50, 92, 82, 12,
        "Real RNA-seq gene-expression data\n(157 samples \u00d7 1,783 genes)\nPreprocessing: quantile normalization,\nlog2 transformation, CoV filtering",
        fill = BOOKEND_FILL, color = BOOKEND_BORDER, textcolor = GREY_TXT)

## ---------- Stage 2: Synthetic generation ----------
add_box(50, 79, 42, 6, "Synthetic-data generation", fill = "white", color = "black", lw = 1.3)

gen_y <- 65; gen_h <- 12
add_box(17, gen_y, 27, gen_h, "ARF\n(tree-based)\nn = 200 synthetic\nsamples", fill = GENERATOR_COLORS["ARF"])
add_box(50, gen_y, 27, gen_h, "CTGAN\n(GAN-based)\nn = 200 synthetic\nsamples", fill = GENERATOR_COLORS["CTGAN"])
add_box(83, gen_y, 27, gen_h, "TVAE\n(VAE-based)\nn = 200 synthetic\nsamples", fill = GENERATOR_COLORS["TVAE"])

## ---------- Stage 3: Evaluation framework outer box ----------
fw_top <- 46; fw_bot <- 23
add_box(50, (fw_top + fw_bot)/2, 94, fw_top - fw_bot, "", fill = FRAMEWORK_FILL,
        color = FRAMEWORK_BORDER, lw = 1.5)

## Non-uniform box widths: "distinguishability" is by far the longest word in
## this row, so it gets more width; short labels (Dilution, etc.) get less.
## This keeps everything on ONE row (matching the previously approved layout)
## while preventing the long word from overflowing its box.
metrics <- c("Predictive\nutility", "Statistical\nfidelity", "Privacy-related\ndistinguishability",
             "Feature-\nimportance\nalignment", "Dilution", "Bootstrap\nstability")
mx <- c(10, 24.5, 44, 63, 76, 88.5)
mw <- c(12, 13, 22, 12, 10, 11)
my <- 31; mh <- 13.5
for (i in seq_along(metrics)) {
  add_box(mx[i], my, mw[i], mh, metrics[i], fill = METRIC_FILL, color = METRIC_BORDER,
          textcolor = "black")
}

## ---------- Stage 4: Comparative analysis ----------
add_box(50, 12.5, 78, 10.5,
        "Comparative analysis across ARF, CTGAN, and TVAE\n(interpreting trade-offs across all six dimensions)",
        fill = BOOKEND_FILL, color = BOOKEND_BORDER, textcolor = GREY_TXT)

rects_df  <- do.call(rbind, rects)
labels_df <- do.call(rbind, labels)

## ---------- arrows ----------
arrows <- data.frame(
  x    = c(50,   38,   50,   62,   17,            50,            83,           50),
  y    = c(86,   76,   76,   76,   gen_y-gen_h/2, gen_y-gen_h/2, gen_y-gen_h/2, fw_bot),
  xend = c(50,   17,   50,   83,   50,            50,            50,           50),
  yend = c(82,   71,   71,   71,   47.5,          47.5,          47.5,         17.8)
)

p <- ggplot() +
  geom_rect(data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rects_df$fill, color = rects_df$color, linewidth = rects_df$lw) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 1.0, color = "black",
               arrow = arrow(length = unit(0.24, "cm"), type = "closed")) +
  geom_text(data = labels_df, aes(x = x, y = y, label = label, color = textcolor),
            family = FONT_FAMILY, fontface = "bold", size = TXT_SIZE, lineheight = 0.95) +
  scale_color_identity() +
  annotate("text", x = 50, y = fw_top - 3.4, label = "Multidimensional evaluation framework",
           family = FONT_FAMILY, fontface = "bold", size = TXT_SIZE, color = FRAMEWORK_BORDER) +
  coord_fixed(xlim = c(0, 100), ylim = c(6, 99.5), expand = FALSE) +
  theme_void()

## Canvas scaled up ~1.15x (matching the 14->16pt font increase) so every box
## gets proportionally more physical room for the same text, without changing
## any box's relative size/position or any text content.
ggsave(file.path(figures_dir, "Fig1_Workflow_overview.png"), p,
       width = 11.5, height = 10.5, dpi = 600, device = ragg::agg_png, bg = "white")
cat("Saved Fig1_Workflow_overview.png\n")
