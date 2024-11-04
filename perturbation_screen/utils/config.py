import os
import re
from typing import Tuple

import yaml
from pydantic import ValidationError
from pydantic.dataclasses import dataclass


@dataclass
class ConfigModel:
    genes_csv: str  # required
    adata_file: str  # required
    group_by: str  # required
    screen_type: str = "DynamoPerturbationScreen"
    perturbation_range: Tuple[int, int] = (-1000, 1000)
    perturbation_steps: int = 21
    perturbation_attempts: int = 2
    working_directory: str = f"{screen_type}_output"
    log_file_template: str = "log-b{rank}.log"
    results_file_template: str = "results-b{rank}.pkl"
    compiled_results_file: str = "results.pkl"
    adata_file_processed: str = ""
    load_processed: bool = False

    def __post_init__(self):
        # Resolve paths
        self.genes_csv = os.path.abspath(self.genes_csv)
        self.adata_file = os.path.abspath(self.adata_file)

        if not self.adata_file_processed:
            filename, ext = os.path.splitext(self.adata_file)
            sanitized_group_by = re.sub(r"[^\w\-_\.]", "_", self.group_by)
            processed_file = "_".join([filename, sanitized_group_by, "processed"]) + ext
            self.adata_file_processed = os.path.abspath(processed_file)


def load_config(config_file: str) -> ConfigModel:
    # Load and validate configuration
    try:
        with open(config_file) as f:
            config_data = yaml.safe_load(f)

        return ConfigModel(**config_data)

    except ValidationError as e:
        print("Configuration error:", e)
        raise
