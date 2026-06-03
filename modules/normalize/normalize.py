import scanpy as sc
import os, sys

in_h5ad = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

print("Step 1: Loading QC data...")
adata = sc.read_h5ad(in_h5ad)
print(f"   Spots: {adata.n_obs} Genes: {adata.n_vars}")

print("Step 2: Library size normalization...")
sc.pp.normalize_total(adata, target_sum=1e4)

print("Step 3: Log transformation...")
sc.pp.log1p(adata)

print("Step 4: Highly Variable Genes selection...")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
print(f" HVGs selected: {adata.var['highly_variable'].sum()}")

print("Step 5: Saving...")
adata.write_h5ad(os.path.join(out_dir, "normalized.h5ad"))
print("Normalization Complete.")
