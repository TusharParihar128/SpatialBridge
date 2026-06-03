import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
import numpy as np

h5_path = sys.argv[1]   # .h5 file
spatial_dir = sys.argv[2]   # spatial/ folder
out_dir = sys.argv[3]   # output folder
os.makedirs(out_dir, exist_ok=True)

print("Ste 1: Loading data...")
adata = sc.read_visium(path=os.path.dirname(h5_path), 
                       count_file = os.path.basename(h5_path))
adata.var_names_make_unique()
print(f" Loaded: {adata.n_obs} spots x {adata.n_vars} genes")

# ── QC Metrics ─────────────────────────────────────────────

print("Step 2: Calculating QC metrics...")
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True,
                           log1p=True)

# ── Adaptive Thresholds (Scran-style) ────────────────────────

print("Step 3: Applying adaptive thresholds...")
def mad(x):
    return np.median(np.abs(x - np.median(x)))
log_counts = np.log1p(adata.obs["total_counts"])
log_genes = np.log1p(adata.obs["n_genes_by_counts"])
mito = adata.obs["pct_counts_mt"]

#Lower Bounds & Upper Bound for Mitochondria
keep = (log_counts >= np.median(log_counts) - 3*mad(log_counts)) & (log_genes >= np.median(log_genes) - 3*mad(log_genes)) & (mito <= np.median(mito) + 3*mad(mito)) & (adata.obs["in_tissue"] == 1)
print(f"    Before: {adata.n_obs}   After: {keep.sum()}     Removed:{(~keep).sum()}")
adata = adata[keep].copy()

# ── Save filtered data ────────────────────────────────────────

print("Step 4: Saving...")
adata.write_h5ad(os.path.join(out_dir, "filtered_qc.h5ad"))

# ── QC Plot ──────────────────────────────────────────────────

print("Step 5: Plotting...")
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(adata.obs["total_counts"], bins=50,
             color="steelblue"); axes[0].set_title("Total Counts per plot")
axes[1].hist(adata.obs["n_genes_by_counts"], bins=50,
             color="seagreen"); axes[1].set_title("Genes per Spot")
axes[2].hist(adata.obs["pct_counts_mt"], bins=50,
             color="tomato"); axes[2].set_title("% Mito Genes")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "qc_plot.png"), dpi=150)
print("QC Complete")

