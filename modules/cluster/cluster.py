import scanpy as sc
import os, sys

in_h5ad = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

print("Step 1: Loading normalized data...")
adata = sc.read_h5ad(in_h5ad)
print(f" Spots: {adata.n_obs} Genes: {adata.n_vars}")

print("Step 2: PCA...")
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

print("Step 3: Neighbourhood graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

print("Step 4: Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
print(f" Clusters found: {adata.obs['leiden'].nunique()}")

print("Step 5: UMAP...")
sc.tl.umap(adata)

print("Step 6: Saving...")
adata.write_h5ad(os.path.join(out_dir, "clustered.h5ad"))

print("Step 7: Plotting...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sc.pl.umap(adata, color='leiden', show=False)
plt.savefig(os.path.join(out_dir, "clusters_umap.png"), dpi=150, bbox_inches='tight')
print("Clustering Complete.")

