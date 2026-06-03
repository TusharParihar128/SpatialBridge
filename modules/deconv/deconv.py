import scanpy as sc
import pandas as pd
import numpy as np
import os, sys

spatial_h5ad = sys.argv[1]
ref_dir      = sys.argv[2]
out_dir      = sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

print("Step 1: Loading spatial data...")
adata = sc.read_h5ad(spatial_h5ad)
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
print(f"  Spots: {adata.n_obs}")

print("Step 2: Loading reference data...")
ref = sc.read_mtx(os.path.join(ref_dir, "count_matrix_sparse.mtx")).T
genes = pd.read_csv(os.path.join(ref_dir, "count_matrix_genes.tsv"), header=None)
barcodes = pd.read_csv(os.path.join(ref_dir, "count_matrix_barcodes.tsv"), header=None)
meta = pd.read_csv(os.path.join(ref_dir, "metadata.csv"), index_col=0)
ref.var_names = genes[0].values
ref.obs_names = barcodes[0].values
ref.obs["celltype_major"] = meta.loc[ref.obs_names, "celltype_major"]
ref.layers["counts"] = ref.X.copy()
print(f"  Reference cells: {ref.n_obs}")
print(f"  Cell types: {ref.obs['celltype_major'].unique().tolist()}")

print("Step 3: Finding common genes...")
common = list(adata.var_names.intersection(ref.var_names))
adata  = adata[:, common].copy()
ref    = ref[:, common].copy()
print(f"  Common genes: {len(common)}")

print("Step 4: Setting up Cell2location reference model...")
import cell2location
from cell2location.models import RegressionModel

sc.pp.filter_genes(ref, min_cells=10)
common2 = list(adata.var_names.intersection(ref.var_names))
adata   = adata[:, common2].copy()
ref     = ref[:, common2].copy()
print(f"  Genes after filtering: {len(common2)}")

RegressionModel.setup_anndata(
    ref,
    labels_key="celltype_major",
    layer="counts"
)

reg_model = RegressionModel(ref)
print("Step 5: Training reference model (200 epochs)...")
reg_model.train(
    max_epochs=200,
    batch_size=2500,
    accelerator="cpu"
)

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
print(f"  Signatures shape: {inf_aver.shape}")

print("Step 7: Setting up Cell2location spatial model...")
cell2location.models.Cell2location.setup_anndata(
    adata,
    layer="counts"
)

c2l_model = cell2location.models.Cell2location(
    adata,
    cell_state_df=inf_aver,
    N_cells_per_location=8,
    detection_alpha=20
)

print("Step 8: Training spatial model (300 epochs)...")
c2l_model.train(
    max_epochs=300,
    batch_size=None,
    accelerator="cpu"
)

print("Step 9: Exporting posterior...")
adata = c2l_model.export_posterior(
    adata,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": adata.n_obs,
        "accelerator": "cpu"
    }
)

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

