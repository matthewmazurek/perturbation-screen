import importlib
import logging
import re
from abc import ABC, abstractmethod
from typing import Any, List, Tuple

import anndata
import numpy as np


class PerturbationScreen(ABC):
    def __init__(self, group_label: str, logger: logging.Logger | None = None):
        self.group_label = group_label
        self.logger = logger

    @abstractmethod
    def load_data(self, adata_file: str) -> anndata.AnnData:
        pass

    def save_data(self, adata: anndata.AnnData, adata_file: str) -> None:
        pass

    @abstractmethod
    def preprocess(self, adata: anndata.AnnData) -> None:
        pass

    @abstractmethod
    def get_group_labels(self, adata: anndata.AnnData) -> List[str]:
        pass

    @abstractmethod
    def get_valid_genes(
        self, adata: anndata.AnnData, gene_list: List[str] | None = None
    ) -> List[str]:
        pass

    @abstractmethod
    def run_experiment(
        self, adata: anndata.AnnData, genes: List[str], expr_vals: List[int]
    ) -> Any:
        pass


def camel_to_snake(name):
    """Convert CamelCase to snake_case (under_score_notation)."""
    return re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()


def import_screen_implementation(implementation_name, base_module="implementations"):
    """Dynamically import the PerturbationScreen class based on the config."""
    module_name = camel_to_snake(implementation_name)
    module_path = f"{base_module}.{module_name}"
    module = importlib.import_module(module_path)
    return getattr(module, implementation_name)


def filter_genes(
    genes: List[str], filter_list: List[str]
) -> Tuple[List[str], List[str]]:
    g, f = np.array(genes), np.array(filter_list)
    mask = np.isin(g, f)
    included, excluded = list(g[mask]), list(g[~mask])
    return included, excluded
