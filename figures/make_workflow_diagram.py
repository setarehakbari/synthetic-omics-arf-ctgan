import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm

plt.rcParams['font.family'] = 'Liberation Serif'
plt.rcParams['font.weight'] = 'bold'

ARF_COLOR   = "#F8766D"
CTGAN_COLOR = "#00BFC4"
TVAE_COLOR  = "#E5A500"
GREY        = "#4D4D4D"
LIGHT_GREY  = "#EAEAEA"

fig, ax = plt.subplots(figsize=(13, 11))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

def box(x, y, w, h, text, facecolor="white", edgecolor="black", fontsize=12,
        textcolor="black", lw=1.6, boxstyle="round,pad=0.02,rounding_size=0.02"):
    b = FancyBboxPatch((x - w/2, y - h/2), w, h,
                        boxstyle=boxstyle,
                        linewidth=lw, edgecolor=edgecolor, facecolor=facecolor, zorder=3)
    ax.add_patch(b)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
             fontweight="bold", color=textcolor, zorder=4, linespacing=1.4)
    return (x, y, w, h)

def varrow(x, y1, y2, lw=2.0, color="black"):
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color,
                                 shrinkA=0, shrinkB=0, mutation_scale=22), zorder=2)

def diag_arrow(x1, y1, x2, y2, lw=1.8, color="black"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color,
                                 shrinkA=2, shrinkB=2, mutation_scale=18), zorder=2)

## ---------- Stage 1: Real data ----------
box(50, 94, 62, 9,
    "Real RNA-seq gene-expression data\n(157 samples \u00d7 1,783 genes)\nPreprocessing: label merge \u2022 filtering \u2022 normalization",
    facecolor=LIGHT_GREY, fontsize=12.5)

varrow(50, 89.5, 83.5)

## ---------- Stage 2: Synthetic generation ----------
box(50, 80, 40, 5.5, "Synthetic-data generation", facecolor="white", edgecolor="black", fontsize=13)

diag_arrow(38, 77.3, 17, 71)
diag_arrow(50, 77.3, 50, 71)
diag_arrow(62, 77.3, 83, 71)

gen_y = 66
gen_h = 11
box(17, gen_y, 26, gen_h, "ARF\n(tree-based)\nn = 200 synthetic\nsamples", facecolor=ARF_COLOR, fontsize=12)
box(50, gen_y, 26, gen_h, "CTGAN\n(GAN-based)\nn = 200 synthetic\nsamples", facecolor=CTGAN_COLOR, fontsize=12)
box(83, gen_y, 26, gen_h, "TVAE\n(VAE-based)\nn = 200 synthetic\nsamples", facecolor=TVAE_COLOR, fontsize=12)

diag_arrow(17, gen_y - gen_h/2, 50, 48.5)
diag_arrow(50, gen_y - gen_h/2, 50, 48.5)
diag_arrow(83, gen_y - gen_h/2, 50, 48.5)

## ---------- Stage 3: Evaluation framework ----------
fw_top, fw_bot = 46.5, 25
box(50, (fw_top+fw_bot)/2, 92, fw_top-fw_bot,
    "", facecolor="#FBFBFB", edgecolor=GREY, lw=1.3, boxstyle="round,pad=0.015,rounding_size=0.02")
ax.text(50, fw_top - 3.2, "Multidimensional evaluation framework",
         ha="center", va="center", fontsize=13.5, fontweight="bold", color=GREY, zorder=4)

metrics = ["Predictive\nutility", "Statistical\nfidelity", "Privacy-related\ndistinguishability",
           "Feature-importance\nalignment", "Dilution", "Bootstrap\nstability"]
mx = [12.5, 28, 43.5, 59, 74.5, 89]
my = 32.5
mw, mh = 14.5, 12.5
for x, label in zip(mx, metrics):
    box(x, my, mw, mh, label, facecolor="white", edgecolor="black", fontsize=10.7)

varrow(50, fw_bot, 20.5, color=GREY)

## ---------- Stage 4: Comparative analysis ----------
box(50, 14, 66, 9.5,
    "Comparative analysis across ARF, CTGAN, and TVAE\n(interpreting trade-offs across all six dimensions)",
    facecolor=LIGHT_GREY, fontsize=12.5)

plt.tight_layout()
plt.savefig("workflow_diagram.png", dpi=600, bbox_inches="tight", facecolor="white")
print("saved")
