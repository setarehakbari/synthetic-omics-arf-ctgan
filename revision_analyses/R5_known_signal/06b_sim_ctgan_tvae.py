"""
06b_sim_ctgan_tvae.py   (known-signal benchmark, step 2)

One CTGAN and one TVAE synthetic data set (200 subjects) per simulated data set,
with exactly the settings of the manuscript:
  CTGAN  50 epochs, batch 256, PAC 1 (50 updates at n = 157)
  TVAE   300 epochs, batch 64, embedding 128, (128,128), L2 1e-5, loss factor 2 (900 updates)
Finished sims are skipped; files are written atomically (restartable).

Usage (from this folder; on a Mac use "caffeinate -i")
  caffeinate -i python3 06b_sim_ctgan_tvae.py
  python3 06b_sim_ctgan_tvae.py --first 1 --last 2 --models CTGAN
  python3 06b_sim_ctgan_tvae.py --quick --last 1          # smoke test (separate folder, NOT usable)
Output: draws/<MODEL>/sim_XX.csv (+ log_XX.csv). Simulated data: no privacy restrictions.
"""
import argparse, math, os, random, time
import numpy as np, pandas as pd, torch
from sdv.metadata import Metadata
from sdv.single_table import CTGANSynthesizer, TVAESynthesizer

CTGAN_CFG = dict(epochs=50, batch_size=256, pac=1)
TVAE_CFG = dict(epochs=300, batch_size=64, embedding_dim=128, compress_dims=(128, 128), decompress_dims=(128, 128),
                l2scale=1e-5, loss_factor=2, enforce_min_max_values=True, enforce_rounding=False)


def set_all_seeds(seed):
    random.seed(seed); np.random.seed(seed % (2 ** 32 - 1)); torch.manual_seed(seed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks-dir", default=os.environ.get("KS_DIR", "."))
    ap.add_argument("--first", type=int, default=1); ap.add_argument("--last", type=int, default=10)
    ap.add_argument("--models", nargs="*", default=["CTGAN", "TVAE"], choices=["CTGAN", "TVAE"])
    ap.add_argument("--n-synth", type=int, default=200)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    seeds = pd.read_csv(os.path.join(args.ks_dir, "seeds.csv"))
    sub = "draws_QUICK" if args.quick else "draws"
    import sdv
    for i in range(args.first, args.last + 1):
        real = pd.read_csv(os.path.join(args.ks_dir, "sims", f"sim_{i:02d}.csv"))
        assert real.shape[0] == 157 and set(real["class"].unique()) == {0, 1} and np.isfinite(real.drop(columns="class").to_numpy()).all()
        metadata = Metadata.detect_from_dataframe(real)
        metadata.update_column("class", sdtype="categorical")
        for name in args.models:
            mdir = os.path.join(args.ks_dir, sub, name); os.makedirs(mdir, exist_ok=True)
            out_csv = os.path.join(mdir, f"sim_{i:02d}.csv")
            if os.path.exists(out_csv):
                print(f"{name} sim {i:02d}: exists, skipping", flush=True); continue
            seed = int(seeds[(seeds.sim == i) & (seeds.role == f"gen_{name}")]["seed"].iloc[0])
            set_all_seeds(seed)
            if name == "CTGAN":
                epochs = 2 if args.quick else CTGAN_CFG["epochs"]
                synth = CTGANSynthesizer(metadata, epochs=epochs, batch_size=CTGAN_CFG["batch_size"], pac=CTGAN_CFG["pac"], verbose=False, enable_gpu=False)
                updates = epochs * max(len(real) // CTGAN_CFG["batch_size"], 1)
            else:
                cfg = dict(TVAE_CFG); cfg["epochs"] = 2 if args.quick else cfg["epochs"]
                synth = TVAESynthesizer(metadata, verbose=False, enable_gpu=False, **cfg)
                updates = cfg["epochs"] * math.ceil(len(real) / cfg["batch_size"])
            print(f"{name} sim {i:02d}: fitting ({updates} updates) ...", flush=True)
            t0 = time.time(); synth.fit(real); fit_s = time.time() - t0
            syn = synth.sample(num_rows=args.n_synth)
            assert set(syn.columns) == set(real.columns)
            syn = syn[list(real.columns)]; syn["class"] = syn["class"].astype(int)
            assert syn.notna().all().all() and set(syn["class"].unique()) <= {0, 1}
            assert np.isfinite(syn.drop(columns="class").to_numpy()).all()
            tmp = out_csv + ".part"; syn.to_csv(tmp, index=False); os.replace(tmp, out_csv)
            pd.DataFrame([dict(sim=i, generator=name, n_synth=len(syn), synth_class0=int((syn["class"] == 0).sum()),
                               synth_class1=int((syn["class"] == 1).sum()), total_updates=updates, fit_seconds=round(fit_s, 1),
                               quick=args.quick, sdv_version=sdv.__version__, torch_version=torch.__version__)]
                         ).to_csv(os.path.join(mdir, f"log_{i:02d}.csv"), index=False)
            print(f"{name} sim {i:02d}: done in {fit_s:.0f} s (class 0/1 = {int((syn['class']==0).sum())}/{int((syn['class']==1).sum())})", flush=True)
    print("CTGAN/TVAE finished.", flush=True)


if __name__ == "__main__":
    main()
