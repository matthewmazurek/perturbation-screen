import hashlib
from pathlib import Path
from typing import Any, List

import anndata
import numpy as np
from interfaces.perturbation_screen import PerturbationScreen
from numpy.typing import NDArray
from utils.logger import log


class DummyPerturbationScreen(PerturbationScreen):
    def load_data(self, adata_file: str) -> anndata.AnnData:
        # Return dummy data instead of actually loading
        log(f"Dummy loading data from {adata_file}...", logger=self.logger)

        # Dummy AnnData with random values
        dummy_adata = anndata.AnnData(np.random.rand(100, 10))

        return dummy_adata

    def save_data(self, adata, adata_file) -> None:
        log(f"Saving dummy preprocessed data to {adata_file}...", logger=self.logger)
        adata.write_h5ad(Path(adata_file))

    def get_group_labels(self, adata: anndata.AnnData) -> List[str]:
        dummy_group_labels = ["group1", "group2", "group3"]
        return dummy_group_labels

    def preprocess(self, adata: anndata.AnnData) -> None:
        # Dummy preprocessing step
        log("Dummy preprocessing step...", logger=None)
        # No actual preprocessing

    def get_valid_genes(
        self, adata: anndata.AnnData, gene_list: List[str] | None = None
    ) -> List[str]:
        # Return dummy gene names
        log("Dummy get_valid_genes call...", logger=None)
        return gene_list or [""]

    def run_experiment(
        self, adata: anndata.AnnData, genes: List[str], expr_vals: List[int]
    ) -> Any:
        # Log the input genes and expression values
        log(
            f"Dummy run_experiment with genes {genes} and expr_vals {expr_vals}...",
            logger=None,
        )

        # Extract the single gene and expression value
        gene = genes[0]
        expr_val = expr_vals[0]

        # Define the number of states
        num_states = len(self.get_group_labels(adata))

        # Initialize the base transition probability matrix
        base_transition_matrix = np.full((num_states, num_states), 1 / num_states)

        # Generate a consistent direction matrix for the gene
        def gene_direction_matrix(gene_name, num_states):
            # Generate a random seed from the gene name
            gene_hash = hashlib.md5(gene_name.encode()).hexdigest()
            seed = int(gene_hash[:8], 16)
            rng = np.random.default_rng(seed)

            # Generate a random direction matrix with values between -1 and 1
            direction_matrix = rng.uniform(-1, 1, size=(num_states, num_states))
            # Ensure adjustments per row sum to zero
            direction_matrix -= direction_matrix.mean(axis=1, keepdims=True)

            return direction_matrix

        direction_matrix = gene_direction_matrix(gene, num_states)

        def adjust_transition_matrix(expr_val, base_matrix, direction_matrix):
            # Compute adjustment factor using tanh to get values between -1 and 1
            # expr_val ∈ [-1000, 1000]
            adjustment_factor = np.tanh(expr_val / 500.0)
            total_adjustment = adjustment_factor * direction_matrix
            adjusted_matrix = base_matrix + total_adjustment

            # Ensure probabilities are valid (non-negative)
            adjusted_matrix = np.clip(adjusted_matrix, 0, None)

            # Normalize the rows to ensure they sum to 1
            adjusted_matrix = adjusted_matrix / adjusted_matrix.sum(
                axis=1, keepdims=True
            )

            return adjusted_matrix

        # Adjust the transition matrix
        transition_matrix = adjust_transition_matrix(
            expr_val, base_transition_matrix, direction_matrix
        )

        return transition_matrix
