import os
from tempfile import NamedTemporaryFile

import pytest
import yaml
from pydantic import ValidationError

from perturbation_screen.utils.config import ConfigModel, load_config


# Test 1: Validation of Required Fields
def test_config_model_missing_required_fields():
    with pytest.raises(ValidationError):
        # Only partial configuration, missing required fields
        ConfigModel(genes_csv="genes.csv", adata_file="data.h5ad")  # type: ignore


# Test 2: Path Resolution


def test_path_resolution():
    config = ConfigModel(
        genes_csv="genes.csv", adata_file="data.h5ad", group_by="group"
    )
    assert os.path.isabs(config.genes_csv)
    assert os.path.isabs(config.adata_file)
    assert os.path.isabs(config.adata_file_processed)  # type: ignore


# Test 3: Sanitization of group_by


def test_group_by_sanitization():
    config = ConfigModel(
        genes_csv="genes.csv", adata_file="data.h5ad", group_by="invalid/group*name"
    )
    # The group_by should be sanitized in adata_file_processed
    assert "invalid_group_name" in config.adata_file_processed  # type: ignore


# Test 4: Default Values


def test_default_values():
    config = ConfigModel(
        genes_csv="genes.csv", adata_file="data.h5ad", group_by="group"
    )
    assert config.screen_type == "DynamoPerturbationScreen"
    assert config.perturbation_range == (-1000, 1000)
    assert config.perturbation_steps == 21
    assert config.perturbation_attempts == 2
    assert config.working_directory == "DynamoPerturbationScreen_output"
    assert config.log_file_template == "log-b{rank}.log"
    assert config.results_file_template == "results-b{rank}.pkl"
    assert config.compiled_results_file == "results.pkl"


# Test 5: Configuration Loading from YAML


def test_load_config():
    config_data = {
        "genes_csv": "genes.csv",
        "adata_file": "data.h5ad",
        "group_by": "group",
    }
    with NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp_file:
        yaml.dump(config_data, tmp_file)
        tmp_file_name = tmp_file.name

    try:
        config = load_config(tmp_file_name)
        assert config.genes_csv == os.path.abspath("genes.csv")
        assert config.adata_file == os.path.abspath("data.h5ad")
        assert config.group_by == "group"
    finally:
        os.remove(tmp_file_name)


# Test 6: Loading Invalid Configuration


def test_load_invalid_config():
    config_data = {
        "genes_csv": "genes.csv",
        # Missing adata_file and group_by
    }
    with NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp_file:
        yaml.dump(config_data, tmp_file)
        tmp_file_name = tmp_file.name

    try:
        with pytest.raises(ValidationError):
            load_config(tmp_file_name)
    finally:
        os.remove(tmp_file_name)
