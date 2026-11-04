"""
05b_draws_ctgan_tvae.py   (repeated-draws pipeline, step 2)

Generates repeated, independent CTGAN / TVAE draws from the FULL real data
(157 subjects x 1,783 genes). One draw = new random seed + new model fit +
new sample of 200 synthetic subjects (the whole procedure is repeated).

Settings (as in the manuscript)
  CTGAN       50 epochs, batch 256, PAC 1
  TVAE        300 epochs, batch 64, embedding 128, (128,128) encoder/decoder,
              L2 = 1e-5, loss factor 2
  CTGAN_long  (optional, slow) CTGAN with as many parameter updates as TVAE
              (n = 157: CTGAN 1 update/epoch, TVAE 3 updates/epoch -> 900 vs 900)

The loop runs draw by draw (draw 1 for every model, then draw 2, ...), so if
you stop early all models have the same number of draws. Finished draws are
skipped; files are written atomically, so the script can be stopped and
restarted at any time.

Usage (run from this folder, with "caffeinate -i" on a Mac)
  caffeinate -i python3 05b_draws_ctgan_tvae.py                    # draws 1-30, CTGAN + TVAE
  python3 05b_draws_ctgan_tvae.py --first 1 --last 5 --models TVAE
  python3 05b_draws_ctgan_tvae.py --first 1 --last 5 --models CTGAN_long   # optional
  python3 05b_draws_ctgan_tvae.py --quick --last 2                  # smoke test (separate folder)

Output: <base>/reviewer_revision/R4_repeated_draws/draws/<MODEL>/draw_XX.csv (+ log_XX.csv)
These files contain synthetic profiles derived from patient data: keep them
OUT of any public repository.
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
    raise ValueError(f"unknown model {name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", default=os.environ.get("RD_BASE_DIR", os.environ.get("SYNTH_EVAL_BASE_DIR", ".")))
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int, default=30)
    ap.add_argument("--models", nargs="*", default=["CTGAN", "TVAE"], choices=["CTGAN", "TVAE", "CTGAN_long"])
    ap.add_argument("--n-synth", type=int, default=200)
    ap.add_argument("--quick", action="store_true", help="2-epoch smoke test; results go to a separate folder and are NOT usable")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    rd_dir = os.path.join(args.base_dir, "reviewer_revision", "R4_repeated_draws")
    seeds = pd.read_csv(os.path.join(rd_dir, "draw_seeds.csv"))
    real = pd.read_csv(os.path.join(args.base_dir, "gene_with_label.csv"))
    assert real.shape == (157, 1784), f"unexpected real data shape {real.shape}"
    assert set(real["class"].unique()) == {0, 1}, "class must be 0/1"
    assert np.isfinite(real.drop(columns="class").to_numpy()).all()
    n_train = len(real)
    metadata = Metadata.detect_from_dataframe(real)
    metadata.update_column("class", sdtype="categorical")      # as in the manuscript scripts

    import sdv
    try:
        import ctgan
        ctgan_version = ctgan.__version__
    except Exception:
        ctgan_version = "NA"
    sub = "draws_QUICK" if args.quick else "draws"

    for d in range(args.first, args.last + 1):
        for name in args.models:
            mdir = os.path.join(rd_dir, sub, name)
            os.makedirs(mdir, exist_ok=True)
            out_csv = os.path.join(mdir, f"draw_{d:02d}.csv")
            if os.path.exists(out_csv) and not args.overwrite:
                print(f"{name} draw {d:02d}: exists, skipping", flush=True)
                continue
            row = seeds[(seeds.generator == name) & (seeds.draw == d) & (seeds.role == "gen")]
            if len(row) != 1:
                raise SystemExit(f"No seed for {name} draw {d}. Run 05_prepare_repeated_draws.R first.")
            seed = int(row["seed"].iloc[0])
            set_all_seeds(seed)
            synth, budget = build_synthesizer(name, metadata, n_train, args.quick)
            print(f"{name} draw {d:02d}: fitting ({budget['epochs']} epochs x {budget['steps_per_epoch']} "
                  f"updates/epoch = {budget['total_updates']} updates) ...", flush=True)
            t0 = time.time()
            synth.fit(real)
            fit_s = time.time() - t0
            syn = synth.sample(num_rows=args.n_synth)

            assert set(syn.columns) == set(real.columns), "synthetic columns differ from the real columns"
            syn = syn[list(real.columns)]
            syn["class"] = syn["class"].astype(int)
            assert syn.notna().all().all() and set(syn["class"].unique()) <= {0, 1}
            assert np.isfinite(syn.drop(columns="class").to_numpy()).all()
            tmp = out_csv + ".part"
            syn.to_csv(tmp, index=False)
            os.replace(tmp, out_csv)                              # atomic: no half-written files
            pd.DataFrame([dict(
                draw=d, generator=name, seed_role="gen", n_real=n_train, n_synth=len(syn),
                synth_class0=int((syn["class"] == 0).sum()), synth_class1=int((syn["class"] == 1).sum()),
                **budget, fit_seconds=round(fit_s, 1), quick=args.quick,
                sdv_version=sdv.__version__, ctgan_version=ctgan_version, torch_version=torch.__version__,
            )]).to_csv(os.path.join(mdir, f"log_{d:02d}.csv"), index=False)
            print(f"{name} draw {d:02d}: done in {fit_s:.0f} s "
                  f"(synthetic class 0/1 = {int((syn['class']==0).sum())}/{int((syn['class']==1).sum())})", flush=True)
    print("CTGAN/TVAE finished.", flush=True)


if __name__ == "__main__":
    main()
