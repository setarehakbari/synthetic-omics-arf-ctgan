# Writes the R session information (R version, package versions) next to the results.
# Run it in the same R installation that was used for the analyses:  Rscript session_info.R
pk <- c("ranger", "arf", "glmnet", "pROC", "rsample", "tidyverse", "npdro", "umap", "patchwork", "cowplot", "ragg")
for (p in pk) suppressWarnings(try(suppressMessages(library(p, character.only = TRUE)), silent = TRUE))
sink("session_info_R.txt"); print(sessionInfo()); sink()
cat("written: session_info_R.txt\n")
# Python (generators): python3 -c "import sdv, ctgan, torch, pandas, numpy; print('sdv', sdv.__version__, 'ctgan', ctgan.__version__, 'torch', torch.__version__, 'pandas', pandas.__version__, 'numpy', numpy.__version__)" > session_info_python.txt
