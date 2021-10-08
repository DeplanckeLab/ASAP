from scopeloompy import *
loom = Loom(file_path="Anndata.loom", title="Test", hierarchy=["","",""])
# from a .h5ad file
loom.from_h5ad(file_path="~/Test/Anndata.h5ad")
# Add embeddings
loom.add_embedding_from_col_attrs(key="X_pca", name="PCA", is_default=True)
loom.add_embedding_from_col_attrs(key="X_umap", name="UMAP", is_default=False)
# Add annotations (discrete variables):
loom.add_annotation_from_col_attrs(key="louvain")
# Add metrics (continuous variables):
loom.add_metric_from_col_attrs(key="n_counts_all")
# Find all the variables that you can add to the loom:
#loom.col_attrs.keys()
# Save to .loom file:
loom.finalize()