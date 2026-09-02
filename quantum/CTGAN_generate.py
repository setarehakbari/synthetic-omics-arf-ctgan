 ============================================================
# CTGAN synthetic gene-expression generation
#
# Input:
#   gene_with_label.csv
#
# Output:
#   synthetic_ctgan_regenerated.csv
#
# This script does NOT overwrite the archived
# synthetic_ctgan.csv used in the manuscript analyses.
# ============================================================

import os
import sys
import time
import random

import numpy as np
import pandas as pd
import torch
import sdv

from sdv.metadata import Metadata
from sdv.single_table import CTGANSynthesizer


# ============================================================
# 1. Reproducibility
# ============================================================

SEED = 42

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)


# ============================================================
# 2. Paths
# ============================================================

base_dir = os.path.expanduser("~/Desktop/My paper")

in_csv = os.path.join(
    base_dir,
    "gene_with_label.csv"
)

out_csv = os.path.join(
    base_dir,
    "synthetic_ctgan_regenerated.csv"
)

meta_json = os.path.join(
    base_dir,
    "ctgan_metadata.json"
)

loss_csv = os.path.join(
    base_dir,
    "ctgan_training_loss.csv"
)


# ============================================================
# 3. Software information
# ============================================================

print("\n========================================")
print("SOFTWARE INFORMATION")
print("========================================")

print("Python version:", sys.version.split()[0])
print("SDV version:", sdv.__version__)
print("Pandas version:", pd.__version__)
print("NumPy version:", np.__version__)
print("PyTorch version:", torch.__version__)
print("Random seed:", SEED)


# ============================================================
# 4. Load real data
# ============================================================

if not os.path.exists(in_csv):
    raise FileNotFoundError(
        f"Input file not found:\n{in_csv}"
    )

df = pd.read_csv(in_csv)

print("\n========================================")
print("REAL DATA")
print("========================================")

print("Input file:", in_csv)
print("Shape:", df.shape)


# ============================================================
# 5. Validate input data
# ============================================================

if "class" not in df.columns:
    raise ValueError(
        "The dataset must contain a column named 'class'."
    )

if df.columns.duplicated().any():
    duplicated = df.columns[
        df.columns.duplicated()
    ].tolist()

    raise ValueError(
        f"Duplicate column names found: {duplicated}"
    )

if df.isna().any().any():
    n_missing = int(
        df.isna().sum().sum()
    )

    raise ValueError(
        f"The dataset contains {n_missing} missing values."
    )


gene_columns = [
    column
    for column in df.columns
    if column != "class"
]


for column in gene_columns:
    df[column] = pd.to_numeric(
        df[column],
        errors="raise"
    )


gene_matrix = df[
    gene_columns
].to_numpy(dtype=float)

if not np.isfinite(
    gene_matrix
).all():

    raise ValueError(
        "The gene-expression matrix contains "
        "Inf or -Inf values."
    )


print("Number of samples:", len(df))
print(
    "Number of gene features:",
    len(gene_columns)
)

print("\nClass counts:")
print(
    df["class"]
    .value_counts()
    .sort_index()
)

print("\nClass proportions:")
print(
    df["class"]
    .value_counts(normalize=True)
    .sort_index()
    .round(3)
)


# ============================================================
# 6. Build metadata
# ============================================================

print("\n========================================")
print("BUILDING METADATA")
print("========================================")


# No primary key exists in this gene-expression dataset.
metadata = Metadata.detect_from_dataframe(
    data=df,
    table_name="gene_expression",
    infer_keys=None
)


# Explicitly define class as categorical.
metadata.update_column(
    column_name="class",
    sdtype="categorical"
)


# Validate metadata before fitting.
metadata.validate()


# Save metadata.
metadata.save_to_json(
    filepath=meta_json,
    mode="overwrite"
)


print("Metadata validated successfully.")
print(
    "Metadata saved to:",
    meta_json
)


# ============================================================
# 7. CTGAN configuration
# ============================================================

EPOCHS = 50
BATCH_SIZE = 256
PAC = 1
SYNTHETIC_ROWS = 200


print("\n========================================")
print("CTGAN CONFIGURATION")
print("========================================")

print("Epochs:", EPOCHS)
print("Batch size:", BATCH_SIZE)
print("PAC:", PAC)
print(
    "Synthetic rows:",
    SYNTHETIC_ROWS
)
print("GPU enabled:", False)
print("Seed:", SEED)


# ============================================================
# 8. Initialize CTGAN
# ============================================================

synthesizer = CTGANSynthesizer(
    metadata=metadata,
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    pac=PAC,
    verbose=True,
    enable_gpu=False
)


# ============================================================
# 9. Fit CTGAN
# ============================================================

print("\n========================================")
print("TRAINING CTGAN")
print("========================================")

start_time = time.time()

synthesizer.fit(df)

elapsed_time = (
    time.time() - start_time
)

print(
    f"\nTraining completed in "
    f"{elapsed_time:.1f} seconds."
)


# ============================================================
# 10. Save training loss
# ============================================================

loss_values = (
    synthesizer.get_loss_values()
)

loss_values.to_csv(
    loss_csv,
    index=False
)

print(
    "Training-loss history saved to:",
    loss_csv
)


# ============================================================
# 11. Generate synthetic data
# ============================================================

print("\n========================================")
print("SAMPLING SYNTHETIC DATA")
print("========================================")

synthetic_data = (
    synthesizer.sample(
        num_rows=SYNTHETIC_ROWS
    )
)

print(
    "Synthetic shape:",
    synthetic_data.shape
)


# ============================================================
# 12. Validate synthetic data
# ============================================================

if list(
    synthetic_data.columns
) != list(df.columns):

    raise RuntimeError(
        "Real and synthetic column names/order "
        "do not match."
    )


if synthetic_data.isna().any().any():

    n_missing = int(
        synthetic_data
        .isna()
        .sum()
        .sum()
    )

    raise RuntimeError(
        f"Synthetic data contain "
        f"{n_missing} missing values."
    )


synthetic_gene_matrix = (
    synthetic_data[
        gene_columns
    ].to_numpy(dtype=float)
)


if not np.isfinite(
    synthetic_gene_matrix
).all():

    raise RuntimeError(
        "Synthetic gene-expression data "
        "contain Inf or -Inf values."
    )


# ============================================================
# 13. Save regenerated synthetic data
# ============================================================

synthetic_data.to_csv(
    out_csv,
    index=False
)

print(
    "\nRegenerated synthetic data saved to:",
    out_csv
)


# ============================================================
# 14. Compare class distributions
# ============================================================

print("\n========================================")
print("CLASS DISTRIBUTION")
print("========================================")


print("\nReal class counts:")
print(
    df["class"]
    .value_counts()
    .sort_index()
)


print("\nSynthetic class counts:")
print(
    synthetic_data["class"]
    .value_counts()
    .sort_index()
)


print("\nReal class proportions:")
print(
    df["class"]
    .value_counts(normalize=True)
    .sort_index()
    .round(3)
)


print("\nSynthetic class proportions:")
print(
    synthetic_data["class"]
    .value_counts(normalize=True)
    .sort_index()
    .round(3)
)


# ============================================================
# 15. Final report
# ============================================================

print("\n========================================")
print("CTGAN RUN COMPLETE")
print("========================================")

print(
    "Real data shape:",
    df.shape
)

print(
    "Synthetic data shape:",
    synthetic_data.shape
)

print(
    "Number of gene features:",
    len(gene_columns)
)

print(
    "Missing values in synthetic data:",
    int(
        synthetic_data
        .isna()
        .sum()
        .sum()
    )
)

print(
    "Output file:",
    out_csv
)

print(
    "\nArchived synthetic_ctgan.csv "
    "was not overwritten."
)

print("\nDone.")
