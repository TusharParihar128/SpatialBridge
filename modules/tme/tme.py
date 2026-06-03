import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
from scipy.sparse import issparse

spatial_h5ad  = sys.argv[1]
deconv_csv    = sys.argv[2]
out_dir       = sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

print("Step 1: Loading data...")
adata = sc.read_h5ad(spatial_h5ad)
prop_df = pd.read_csv(deconv_csv, index_col=0)
common = adata.obs_names.intersection(prop_df.index)
adata = adata[common].copy()
prop_df = prop_df.loc[common]
for col in prop_df.columns:
    adata.obs[col] = prop_df[col].values
print(f"  Spots: {adata.n_obs}, Cell types: {len(prop_df.columns)}")

print("Step 2: Dominant cell type per spot...")
adata.obs['dominant_cell_type'] = pd.Categorical(prop_df.idxmax(axis=1).values)
print(adata.obs['dominant_cell_type'].value_counts())

print("Step 3: Immune infiltration score...")
immune_cols = [c for c in prop_df.columns if c in ['T-cells','B-cells','Myeloid','Plasmablasts']]
tumor_cols  = [c for c in prop_df.columns if c in ['Cancer Epithelial']]
adata.obs['immune_score'] = prop_df[immune_cols].sum(axis=1).values
adata.obs['tumor_score']  = prop_df[tumor_cols].sum(axis=1).values
print(f"  Mean immune score: {adata.obs['immune_score'].mean():.3f}")
print(f"  Mean tumor score:  {adata.obs['tumor_score'].mean():.3f}")

print("Step 4: Spatial co-localization...")
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
sq.gr.nhood_enrichment(adata, cluster_key="dominant_cell_type")
print("  Neighborhood enrichment done.")

print("Step 5: Moran's I - spatial autocorrelation...")
sq.gr.spatial_autocorr(adata, mode="moran", genes=prop_df.columns.tolist(), use_raw=False, layer=None, attr="obs")
print(adata.uns["moranI"].head())

print("Step 6: Saving plots...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
axes[0].scatter(range(adata.n_obs), adata.obs['immune_score'], s=1, c='steelblue')
axes[0].set_title("Immune Score per Spot")
axes[1].scatter(range(adata.n_obs), adata.obs['tumor_score'], s=1, c='tomato')
axes[1].set_title("Tumor Score per Spot")
axes[2].bar(adata.obs['dominant_cell_type'].value_counts().index,
            adata.obs['dominant_cell_type'].value_counts().values, color='seagreen')
axes[2].set_title("Dominant Cell Types")
axes[2].tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "tme_plots.png"), dpi=150, bbox_inches='tight')

print("Step 7: Saving results...")
adata.obs.to_csv(os.path.join(out_dir, "tme_results.csv"))
adata.uns["moranI"].to_csv(os.path.join(out_dir, "moransI.csv"))
print("TME Analysis Complete.")
