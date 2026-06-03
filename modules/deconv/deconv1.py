import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import torch

# Arguments validation
if len(sys.argv) < 4:
    print("Usage: python deconv.py <spatial_data_path> <ref_dir> <out_dir>")
    sys.exit(1)

spatial_path = sys.argv[1]
ref_dir      = sys.argv[2]
out_dir      = sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

# ==========================================
# Step 1: Loading Spatial Data
# ==========================================
print("Step 1: Loading spatial data...")

if spatial_path.endswith('.h5'):
    print(f"  Detected 10X H5 file format. Loading from: {spatial_path}")
    adata = sc.read_10x_h5(spatial_path)
else:
    print(f"  Detected H5AD format. Loading from: {spatial_path}")
    adata = sc.read_h5ad(spatial_path)

# Ensure unique gene names right after loading
adata.var_names_make_unique()

# --- ABSOLUTE SAFETY NET: FILTER COMPLETELY EMPTY SPOTS ---
# Sirf un spots ko rakho jahan minimum 10 genes aur minimum 100 counts hon (Removes background noise)
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_cells(adata, min_counts=100)

# Cast data matrix strictly to integer (GammaPoisson support requirement)
if hasattr(adata.X, "data"):
    adata.X.data = np.ceil(adata.X.data).astype(np.int64)
else:
    adata.X = np.ceil(adata.X).astype(np.int64)

# Cell2location batch effects tracker requirement
adata.obs['sample'] = 'sample_1'
adata.layers["counts"] = adata.X.copy()
print(f"  Cleaned Spots for training: {adata.n_obs}")

# ==========================================
# Step 2: Loading Reference Data
# ==========================================
print("Step 2: Loading reference data...")
ref = sc.read_mtx(os.path.join(ref_dir, "count_matrix_sparse.mtx")).T
genes = pd.read_csv(os.path.join(ref_dir, "count_matrix_genes.tsv"), header=None)
barcodes = pd.read_csv(os.path.join(ref_dir, "count_matrix_barcodes.tsv"), header=None)
meta = pd.read_csv(os.path.join(ref_dir, "metadata.csv"), index_col=0)

ref.var_names = genes[0].values
ref.obs_names = barcodes[0].values
ref.obs["celltype_major"] = meta.loc[ref.obs_names, "celltype_major"]

ref.var_names_make_unique()

# Clean reference cells as well to remove absolute zeros
sc.pp.filter_cells(ref, min_genes=10)
sc.pp.filter_cells(ref, min_counts=100)

if hasattr(ref.X, "data"):
    ref.X.data = np.ceil(ref.X.data).astype(np.int64)
else:
    ref.X = np.ceil(ref.X).astype(np.int64)

ref.layers["counts"] = ref.X.copy()
print(f"  Cleaned Reference cells: {ref.n_obs}")
print(f"  Cell types: {ref.obs['celltype_major'].unique().tolist()}")

# ==========================================
# Step 3: Finding Common Genes
# ==========================================
print("Step 3: Finding common genes...")
common = list(adata.var_names.intersection(ref.var_names))
adata  = adata[:, common].copy()
ref    = ref[:, common].copy()
print(f"  Common genes matching matrix: {len(common)}")

# ==========================================
# Step 4: Setting up Cell2location Reference Model
# ==========================================
print("Step 4: Setting up Cell2location reference model...")
import cell2location
from cell2location.models import RegressionModel

sc.pp.filter_genes(ref, min_cells=15)
common2 = list(adata.var_names.intersection(ref.var_names))
adata   = adata[:, common2].copy()
ref     = ref[:, common2].copy()
print(f"  Genes after stringent filtering: {len(common2)}")

RegressionModel.setup_anndata(
    ref,
    labels_key="celltype_major",
    layer="counts"
)

reg_model = RegressionModel(ref)

# ==========================================
# Step 5: Training Reference Model
# ==========================================
print("Step 5: Training reference model (200 epochs)...")
reg_model.train(
    max_epochs=200,
    batch_size=2500,
    accelerator="cpu"
)

# ==========================================
# Step 6: Extracting Cell Type Signatures
# ==========================================
print("Step 6: Extracting cell type signatures...")
ref = reg_model.export_posterior(
    ref,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": 2500,
        "accelerator": "cpu"
    }
)
if "means_per_cluster_mu_fg" in ref.varm:
    inf_aver = pd.DataFrame(
        ref.varm["means_per_cluster_mu_fg"],
        index=ref.var_names,
        columns=ref.uns["mod"]["factor_names"]
    )
else:
    inf_aver = pd.DataFrame(
        ref.obsm["means_per_cluster_mu_fg"],
        index=ref.var_names,
        columns=ref.uns["mod"]["factor_names"]
    )

# --- CRITICAL STABILITY SHIELD: CLIP ABSOLUTE ZEROS FROM SIGNATURES ---
# Agar koi gene kisi cell type mein bilkul express nahi hai, toh use 1e-5 set karo taaki rate check crash na ho
inf_aver = inf_aver.clip(lower=1e-5)
print(f"  Signatures shape: {inf_aver.shape}")

# ==========================================
# Step 7: Setting up Cell2location Spatial Model
# ==========================================
print("Step 7: Setting up Cell2location spatial model...")
cell2location.models.Cell2location.setup_anndata(
    adata,
    batch_key="sample",
    layer="counts"
)

c2l_model = cell2location.models.Cell2location(
    adata,
    cell_state_df=inf_aver,
    N_cells_per_location=8,
    detection_alpha=20
)

# ==========================================
# Step 8: Training Spatial Model
# ==========================================
print("Step 8: Training spatial model (300 epochs)...")
c2l_model.train(
    max_epochs=300,
    batch_size=None,
    accelerator="cpu"
)

# ==========================================
# Step 9: Exporting Posterior
# ==========================================
print("Step 9: Exporting posterior...")
adata = c2l_model.export_posterior(
    adata,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": adata.n_obs,
        "accelerator": "cpu"
    }
)

# ==========================================
# Step 10: Saving Results
# ==========================================
print("Step 10: Saving results...")
adata.write_h5ad(os.path.join(out_dir, "deconvolved.h5ad"))

if "q05_cell_abundance_w_sf" in adata.obsm:
    prop_df = pd.DataFrame(
        adata.obsm["q05_cell_abundance_w_sf"],
        index=adata.obs_names,
        columns=adata.uns["mod"]["factor_names"]
    )
    prop_df.to_csv(os.path.join(out_dir, "cell_proportions.csv"))
    print("\nAverage cell abundances:")
    print(prop_df.mean().sort_values(ascending=False).round(3))

print("\nDeconvolution Complete.")
