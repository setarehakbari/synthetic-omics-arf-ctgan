"""
02b_generate_ctgan_tvae.py   (stage 2b of 3)

Fits CTGAN and TVAE (SDV) on the TRAINING subjects of each split only
and generates 200 synthetic subjects per model.

Three synthesizers are run per split:
  CTGAN       settings of the manuscript (50 epochs, batch 256, PAC = 1)
  CTGAN_long  same CTGAN settings, but with as many generator updates
              as the manuscript TVAE run receives (see below)
  TVAE        settings of the manuscript (300 epochs, batch 64,
              embedding 128, (128,128) encoder/decoder, L2 = 1e-5,
              loss factor 2)

Training budget bookkeeping
  Reviewer 4 observed that CTGAN received far fewer parameter updates
  than TVAE. For n training subjects:
      CTGAN updates per epoch = max(n // batch_size, 1)
      TVAE  updates per epoch = ceil(n / batch_size)
  The log written for every fit records epochs, batch size, updates
  per epoch, total updates and wall-clock training time.
  NOTE: a CTGAN "update" (one discriminator step plus one generator
  step) and a TVAE "update" (one joint optimiser step) are not
  computationally equivalent; matched counts are a reference point,
  not proof of equal training quality.

Usage
  python 02b_generate_ctgan_tvae.py                       # all splits, all models
  python 02b_generate_ctgan_tvae.py --splits 1 2          # selected splits
  python 02b_generate_ctgan_tvae.py --models CTGAN TVAE   # selected models
  python 02b_generate_ctgan_tvae.py --quick               # 2-epoch smoke test only
Existing outputs are kept unless --overwrite is given. Several terminals
can work on different splits at the same time.

Output per split and model:
  splits/split_XX/synthetic_<MODEL>.csv
  splits/split_XX/log_<MODEL>.csv
"""
import argparse
import math
import os
import random
import time

import numpy as np
import pandas as pd
import torch
from sdv.metadata import Metadata
from sdv.single_table import CTGANSynthesizer, TVAESynthesizer

CTGAN_CFG = dict(epochs=50, batch_size=256, pac=1)
TVAE_CFG = dict(epochs=300, batch_size=64, embedding_dim=128,
                compress_dims=(128, 128), decompress_dims=(128, 128),
                l2scale=1e-5, loss_factor=2,
                enforce_min_max_values=True, enforce_rounding=False)


def ctgan_steps_per_epoch(n, batch_size):
    return max(n // batch_size, 1)


def tvae_steps_per_epoch(n, batch_size):
    return math.ceil(n / batch_size)


def set_all_seeds(seed):
    random.seed(seed)
    np.random.seed(seed % (2 ** 32 - 1))
    torch.manual_seed(seed)


def build_synthesizer(name, metadata, n_train, quick):
    """Return (synthesizer, budget dict)."""
    tvae_total = TVAE_CFG["epochs"] * tvae_steps_per_epoch(n_train, TVAE_CFG["batch_size"])
    if name in ("CTGAN", "CTGAN_long"):
        spe = ctgan_steps_per_epoch(n_train, CTGAN_CFG["batch_size"])
        epochs = CTGAN_CFG["epochs"] if name == "CTGAN" else math.ceil(tvae_total / spe)
        if quick:
            epochs = 2 if name == "CTGAN" else 3
        synth = CTGANSynthesizer(metadata, epochs=epochs, batch_size=CTGAN_CFG["batch_size"],
                                 pac=CTGAN_CFG["pac"], verbose=False, enable_gpu=False)
        return synth, dict(epochs=epochs, batch_size=CTGAN_CFG["batch_size"], steps_per_epoch=spe,
                           total_updates=epochs * spe)
    if name == "TVAE":
        spe = tvae_steps_per_epoch(n_train, TVAE_CFG["batch_size"])
        cfg = dict(TVAE_CFG)
        if quick:
            cfg["epochs"] = 2
        synth = TVAESynthesizer(metadata, verbose=False, enable_gpu=False, **cfg)
        return synth, dict(epochs=cfg["epochs"], batch_size=cfg["batch_size"], steps_per_epoch=spe,
                           total_updates=cfg["epochs"] * spe)
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", default=os.environ.get("HOLDOUT_BASE_DIR",
                                                          os.environ.get("SYNTH_EVAL_BASE_DIR", ".")))
    ap.add_argument("--splits", nargs="*", type=int, default=None)
    ap.add_argument("--models", nargs="*", default=["CTGAN", "CTGAN_long", "TVAE"])
    ap.add_argument("--n-synth", type=int, default=200)
    ap.add_argument("--quick", action="store_true", help="2-epoch smoke test; outputs are NOT usable results")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    out_dir = os.path.join(args.base_dir, "reviewer_revision", "R3_major2_holdout")
    seeds = pd.read_csv(os.path.join(out_dir, "seeds.csv"))
    all_splits = sorted(int(d.split("_")[1]) for d in os.listdir(os.path.join(out_dir, "splits"))
                        if d.startswith("split_"))
    splits = args.splits if args.splits else all_splits

    import sdv
    try:
        import ctgan
        ctgan_version = ctgan.__version__
    except Exception:
        ctgan_version = "NA"

    for i in splits:
        d = os.path.join(out_dir, "splits", f"split_{i:02d}")
        train = pd.read_csv(os.path.join(d, "train.csv"))
        n_train = len(train)
        metadata = Metadata.detect_from_dataframe(train)
        metadata.update_column("class", sdtype="categorical")   # as in the manuscript scripts

        for name in args.models:
            out_csv = os.path.join(d, f"synthetic_{name}.csv")
            if os.path.exists(out_csv) and not args.overwrite:
                print(f"split {i:02d} {name}: output exists, skipping", flush=True)
                continue
            seed = int(seeds.loc[(seeds.split == i) & (seeds.role == f"gen_{name}"), "seed"].iloc[0])
            set_all_seeds(seed)
            synth, budget = build_synthesizer(name, metadata, n_train, args.quick)
            print(f"split {i:02d} {name}: fitting ({budget['epochs']} epochs x {budget['steps_per_epoch']} "
                  f"updates/epoch = {budget['total_updates']} updates) ...", flush=True)
            t0 = time.time()
            synth.fit(train)
            fit_s = time.time() - t0
            syn = synth.sample(num_rows=args.n_synth)

            assert set(syn.columns) == set(train.columns), "synthetic columns differ from training columns"
            syn = syn[list(train.columns)]
            syn["class"] = syn["class"].astype(int)
            assert syn.notna().all().all() and set(syn["class"].unique()) <= {0, 1}
            suffix = ".QUICK" if args.quick else ""
            syn.to_csv(out_csv + suffix, index=False)
            pd.DataFrame([dict(
                split=i, generator=name, n_train=n_train, n_synth=len(syn),
                synth_class0=int((syn["class"] == 0).sum()), synth_class1=int((syn["class"] == 1).sum()),
                **budget, fit_seconds=round(fit_s, 1), seed_role=f"gen_{name}", quick=args.quick,
                sdv_version=sdv.__version__, ctgan_version=ctgan_version, torch_version=torch.__version__,
            )]).to_csv(os.path.join(d, f"log_{name}.csv" + suffix), index=False)
            print(f"split {i:02d} {name}: done in {fit_s:.0f} s "
                  f"(synthetic class 0/1 = {int((syn['class']==0).sum())}/{int((syn['class']==1).sum())})", flush=True)


if __name__ == "__main__":
    main()
