############################################################
# TVAE_generate.py
#
# TVAE synthetic gene-expression generation
# Additional generator for reviewer revision
#
# Input:
#   gene_with_label.csv
#
# Output:
#   synthetic_tvae.csv
#
# The script:
#   - validates the real dataset
#   - treats "class" as categorical
#   - trains TVAE on all gene-expression features
#   - generates 200 synthetic samples
#   - preserves column names/order
#   - exports training loss, metadata, and configuration
############################################################

from pathlib import Path
import os
import sys
import json
import time
import random
import platform
import importlib.metadata

import numpy as np
import pandas as pd
import torch

from sdv.metadata import Metadata
from sdv.single_table import TVAESynthesizer


# ==========================================================
# 1. Configuration
# ==========================================================

SEED = 42

BASE_DIR = Path.home() / "Desktop" / "My paper"

INPUT_PATH = BASE_DIR / "gene_with_label.csv"
OUTPUT_PATH = BASE_DIR / "synthetic_tvae.csv"

AUX_DIR = BASE_DIR / "TVAE_generation"
AUX_DIR.mkdir(parents=True, exist_ok=True)

METADATA_PATH = AUX_DIR / "tvae_metadata.json"
LOSS_PATH = AUX_DIR / "tvae_training_loss.csv"
CONFIG_PATH = AUX_DIR / "tvae_run_config.json"

# Expected dataset structure
EXPECTED_GENE_COUNT = 1783

# Use the same synthetic sample count as ARF and CTGAN
SYNTHETIC_ROWS = 200

# TVAE configuration
EPOCHS = 300
BATCH_SIZE = 64

EMBEDDING_DIM = 128
COMPRESS_DIMS = (128, 128)
DECOMPRESS_DIMS = (128, 128)

L2SCALE = 1e-5
LOSS_FACTOR = 2

# CPU is used for consistency across systems.
ENABLE_GPU = False


# ==========================================================
# 2. Random seeds
# ==========================================================

os.environ["PYTHONHASHSEED"] = str(SEED)

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# Ask PyTorch to prefer deterministic operations where possible.
try:
    torch.use_deterministic_algorithms(True, warn_only=True)
except Exception:
    pass


# ==========================================================
# 3. Software versions
# ==========================================================

def package_version(package_name):
    try:
        return importlib.metadata.version(package_name)
    except Exception:
        return "unknown"


software_versions = {
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "pandas": pd.__version__,
    "numpy": np.__version__,
    "torch": torch.__version__,
    "sdv": package_version("sdv"),
}


print("\n====================================================")
print("TVAE synthetic-data generation")
print("====================================================")

print("\nSoftware versions:")
for key, value in software_versions.items():
    print(f"  {key}: {value}")


# ==========================================================
# 4. Load and validate real data
# ==========================================================

if not INPUT_PATH.exists():
    raise FileNotFoundError(
        f"Real data file was not found:\n{INPUT_PATH}"
    )


df = pd.read_csv(INPUT_PATH)

print("\nReal dataset loaded.")
print("Shape:", df.shape)


if "class" not in df.columns:
    raise ValueError(
        "The input dataset must contain a column named 'class'."
    )


gene_columns = [col for col in df.columns if col != "class"]

print("Number of gene features:", len(gene_columns))
print("\nReal class distribution:")
print(df["class"].value_counts(dropna=False).sort_index())


# ----------------------------------------------------------
# Check expected dimensions
# ----------------------------------------------------------

if len(gene_columns) != EXPECTED_GENE_COUNT:
    raise ValueError(
        f"Expected {EXPECTED_GENE_COUNT} gene features, "
        f"but found {len(gene_columns)}."
    )


# ----------------------------------------------------------
# Convert all genes to numeric
# ----------------------------------------------------------

for col in gene_columns:
    df[col] = pd.to_numeric(df[col], errors="raise")


# ----------------------------------------------------------
# Standardize class representation
# ----------------------------------------------------------

class_as_string = df["class"].astype(str).str.strip()

valid_classes = {"0", "1"}

if not set(class_as_string.unique()).issubset(valid_classes):
    raise ValueError(
        "The 'class' column must contain only 0 and 1."
    )

df["class"] = class_as_string.astype(int)


# ----------------------------------------------------------
# Check missing/infinite values
# ----------------------------------------------------------

missing_count = int(df.isna().sum().sum())

if missing_count != 0:
    raise ValueError(
        f"Input dataset contains {missing_count} missing values."
    )


gene_array = df[gene_columns].to_numpy(dtype=float)

if not np.isfinite(gene_array).all():
    raise ValueError(
        "Input gene-expression data contain non-finite values."
    )


duplicated_rows = int(df.duplicated().sum())

print("\nInput validation:")
print("  Rows:", len(df))
print("  Gene features:", len(gene_columns))
print("  Missing values:", missing_count)
print("  Duplicated rows:", duplicated_rows)


# ==========================================================
# 5. Construct SDV metadata
# ==========================================================

metadata = Metadata.detect_from_dataframe(
    data=df,
    table_name="gene_expression",
    infer_keys=None
)

# The diagnosis variable must be categorical.
metadata.update_column(
    column_name="class",
    sdtype="categorical"
)

metadata.validate()

metadata.save_to_json(
    filepath=str(METADATA_PATH),
    mode="overwrite"
)

print("\nMetadata validated and saved:")
print(METADATA_PATH)


# ==========================================================
# 6. Initialize TVAE
# ==========================================================

print("\nTVAE configuration:")
print("  Seed:", SEED)
print("  Epochs:", EPOCHS)
print("  Batch size:", BATCH_SIZE)
print("  Embedding dimension:", EMBEDDING_DIM)
print("  Encoder dimensions:", COMPRESS_DIMS)
print("  Decoder dimensions:", DECOMPRESS_DIMS)
print("  L2 scale:", L2SCALE)
print("  Loss factor:", LOSS_FACTOR)
print("  GPU enabled:", ENABLE_GPU)
print("  Synthetic rows:", SYNTHETIC_ROWS)


synthesizer = TVAESynthesizer(
    metadata=metadata,

    # Gene-expression values are continuous.
    enforce_min_max_values=True,
    enforce_rounding=False,

    embedding_dim=EMBEDDING_DIM,
    compress_dims=COMPRESS_DIMS,
    decompress_dims=DECOMPRESS_DIMS,

    l2scale=L2SCALE,
    batch_size=BATCH_SIZE,
    epochs=EPOCHS,
    loss_factor=LOSS_FACTOR,

    verbose=True,
    enable_gpu=ENABLE_GPU
)


# ==========================================================
# 7. Train TVAE
# ==========================================================

print("\n====================================================")
print("Training TVAE...")
print("====================================================\n")

start_time = time.time()

synthesizer.fit(df)

elapsed_seconds = time.time() - start_time

print("\nTVAE training completed.")
print(
    f"Training time: {elapsed_seconds / 60:.2f} minutes"
)


# ==========================================================
# 8. Save training loss
# ==========================================================

try:
    loss_values = synthesizer.get_loss_values()

    loss_values.to_csv(
        LOSS_PATH,
        index=False
    )

    print("\nTraining loss saved:")
    print(LOSS_PATH)

except Exception as exc:
    print("\nWarning: training loss could not be exported.")
    print(exc)


# ==========================================================
# 9. Generate synthetic data
# ==========================================================

print("\nGenerating synthetic samples...")

# Reset sampling state before generating the manuscript dataset.
try:
    synthesizer.reset_sampling()
except Exception:
    pass


synthetic = synthesizer.sample(
    num_rows=SYNTHETIC_ROWS
)


# ==========================================================
# 10. Validate generated data
# ==========================================================

# Column names must match the real dataset.
if set(synthetic.columns) != set(df.columns):
    missing_columns = set(df.columns) - set(synthetic.columns)
    extra_columns = set(synthetic.columns) - set(df.columns)

    raise ValueError(
        "Synthetic dataset columns do not match real data.\n"
        f"Missing columns: {missing_columns}\n"
        f"Extra columns: {extra_columns}"
    )


# Preserve original column order.
synthetic = synthetic[df.columns]


# Convert gene columns to numeric.
for col in gene_columns:
    synthetic[col] = pd.to_numeric(
        synthetic[col],
        errors="raise"
    )


# Normalize class representation to integer 0/1.
synthetic_class = (
    synthetic["class"]
    .astype(str)
    .str.strip()
)

if not set(synthetic_class.unique()).issubset(valid_classes):
    raise ValueError(
        "TVAE generated unexpected values in the class column: "
        f"{sorted(synthetic_class.unique())}"
    )

synthetic["class"] = synthetic_class.astype(int)


# Check dimensions.
expected_columns = EXPECTED_GENE_COUNT + 1

if synthetic.shape != (SYNTHETIC_ROWS, expected_columns):
    raise ValueError(
        "Unexpected synthetic-data dimensions. "
        f"Expected ({SYNTHETIC_ROWS}, {expected_columns}) "
        f"but obtained {synthetic.shape}."
    )


# Check missing values.
synthetic_missing = int(
    synthetic.isna().sum().sum()
)

if synthetic_missing != 0:
    raise ValueError(
        f"Synthetic data contain "
        f"{synthetic_missing} missing values."
    )


# Check finite gene-expression values.
synthetic_gene_array = synthetic[
    gene_columns
].to_numpy(dtype=float)

if not np.isfinite(synthetic_gene_array).all():
    raise ValueError(
        "Synthetic gene-expression data contain "
        "non-finite values."
    )


# ==========================================================
# 11. Save synthetic dataset
# ==========================================================

synthetic.to_csv(
    OUTPUT_PATH,
    index=False
)


print("\n====================================================")
print("Synthetic TVAE dataset successfully generated")
print("====================================================")

print("\nSynthetic shape:")
print(synthetic.shape)

print("\nSynthetic class distribution:")
print(
    synthetic["class"]
    .value_counts()
    .sort_index()
)

print("\nMissing values:")
print(synthetic_missing)

print("\nDuplicated synthetic rows:")
print(int(synthetic.duplicated().sum()))

print("\nSynthetic dataset saved to:")
print(OUTPUT_PATH)


# ==========================================================
# 12. Save run configuration
# ==========================================================

run_config = {
    "input_file": str(INPUT_PATH),
    "output_file": str(OUTPUT_PATH),

    "random_seed": SEED,

    "real_rows": int(df.shape[0]),
    "gene_features": int(len(gene_columns)),

    "real_class_distribution": {
        str(k): int(v)
        for k, v in
        df["class"].value_counts().sort_index().items()
    },

    "synthetic_rows": SYNTHETIC_ROWS,

    "synthetic_class_distribution": {
        str(k): int(v)
        for k, v in
        synthetic["class"].value_counts().sort_index().items()
    },

    "model": "SDV TVAESynthesizer",

    "epochs": EPOCHS,
    "batch_size": BATCH_SIZE,

    "embedding_dim": EMBEDDING_DIM,
    "compress_dims": list(COMPRESS_DIMS),
    "decompress_dims": list(DECOMPRESS_DIMS),

    "l2scale": L2SCALE,
    "loss_factor": LOSS_FACTOR,

    "enforce_min_max_values": True,
    "enforce_rounding": False,
    "enable_gpu": ENABLE_GPU,

    "training_time_seconds": elapsed_seconds,

    "software_versions": software_versions,
}


with open(
    CONFIG_PATH,
    "w",
    encoding="utf-8"
) as f:

    json.dump(
        run_config,
        f,
        indent=2
    )


print("\nRun configuration saved:")
print(CONFIG_PATH)

print("\nDone.")
