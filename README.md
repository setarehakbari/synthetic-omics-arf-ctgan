# synthetic-omics-arf-ctgan
# Synthetic-omics-arf-ctgan

This repository contains code for evaluating synthetic gene-expression data in major depressive disorder (MDD), with a comparison between Adversarial Random Forest and Conditional Tabular GAN.

## Repository structure
- `R/`: R scripts for figures and tables used in the manuscript
  - `R/Table2_rank_metrics.R`: computes Spearman correlation, Kendall correlation, Rank-Biased Overlap (RBO), Top-50 Jaccard similarity, mean and median absolute rank shifts, and Top-50 dropout for Real vs. Real+ARF and Real vs. Real+CTGAN comparisons.
- `quantum/`: Python scripts for the quantum-adversarial extension
- `results/`: generated figures and tables used in the manuscript

## Data availability
The real gene-expression dataset used in this study is not publicly available due to data-use restrictions. The repository includes code to reproduce the analyses for researchers with authorized access to the data.

