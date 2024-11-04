from pathlib import Path
from typing import Any, List

import anndata
import dynamo as dyn
from interfaces.perturbation_screen import PerturbationScreen
from numpy.typing import NDArray
from utils.logger import log


class DynamoPerturbationScreen(PerturbationScreen):
    def load_data(self, adata_file: str) -> anndata.AnnData:
        log(f"Loading data from {adata_file}...", logger=self.logger)
        adata = dyn.read(adata_file)
        return adata

    def save_data(self, adata, adata_file) -> None:
        log(f"Saving preprocessed data to {adata_file}...", logger=self.logger)
        adata.write_h5ad(Path(adata_file))

    def preprocess(self, adata: anndata.AnnData) -> None:
        log("Preprocessing data...", logger=self.logger)
        pp = dyn.pp.Preprocessor()
        pp.preprocess_adata(adata, recipe="monocle")
        dyn.tl.dynamics(adata)
        dyn.tl.reduceDimension(adata)
        dyn.tl.cell_velocities(adata, basis="pca")
        dyn.vf.VectorField(adata, basis="pca", M=100)
        dyn.tl.cell_velocities(adata, basis="umap")
        dyn.vf.VectorField(adata, basis="umap", M=100)
        dyn.pd.state_graph(adata, group=self.group_label, basis="umap")

    def get_group_labels(self, adata: anndata.AnnData) -> List[str]:
        return list(set(adata.obs[self.group_label]))

    def get_valid_genes(
        self, adata: anndata.AnnData, gene_list: List[str] | None = None
    ) -> List[str]:
        pca_genes = adata.var_names[adata.var.use_for_pca]
        return list(pca_genes)

    def run_experiment(
        self, adata: anndata.AnnData, genes: List[str], expr_vals: List[int]
    ) -> Any:
        dyn.pd.perturbation(adata, genes, expr_vals, emb_basis="umap")
        dyn.vf.VectorField(adata, basis="umap_perturbation")
        dyn.pd.state_graph(
            adata,
            group=self.group_label,
            method="vf",
            basis="umap_perturbation",
            transition_mat_key="perturbation_transition_matrix",
        )

        return adata.uns[f"{self.group_label}_graph"]["group_graph"].copy()
