from unittest import mock

import pytest
import yaml

from perturbation_screen.interfaces.perturbation_screen import (
    import_screen_implementation,
)
from perturbation_screen.main import run, save_results


@pytest.fixture
def mock_makedirs():
    with mock.patch("os.makedirs") as makedirs_mock:
        yield makedirs_mock


@pytest.fixture
def mock_open():
    with mock.patch("builtins.open", mock.mock_open()) as open_mock:
        yield open_mock


@pytest.fixture
def mock_pickle_dump():
    with mock.patch("pickle.dump") as pickle_mock:
        yield pickle_mock


@pytest.fixture
def mock_comm():
    # Mocking MPI communication calls
    mock_comm = mock.Mock()
    mock_comm.Get_rank.return_value = 0
    mock_comm.Get_size.return_value = 4
    yield mock_comm


# Test: save_results


def test_save_results():
    results = {"key": "value"}
    path = "./results.pkl"

    # Patch both open and pickle.dump
    with mock.patch("builtins.open", mock.mock_open()) as mock_open:
        with mock.patch("pickle.dump") as mock_pickle_dump:
            save_results(results, path)

            # Ensure the file was opened in write-binary mode
            mock_open.assert_called_once_with(path, "wb")

            # Get the file handle that would be returned by open()
            file_handle = mock_open()

            # Ensure pickle.dump was called to save the data
            mock_pickle_dump.assert_called_once_with(results, file_handle)


def test_run(mock_open, mock_pickle_dump, mock_comm):
    screen = mock.Mock()  # Mock the PerturbationScreen
    adata = mock.Mock()  # Mock the AnnData object
    genes_list = ["gene1", "gene2"]
    expr_vals = [0.5, 1.0]
    results_file = "./results.pkl"
    rank = 0

    # Simulate starting a new run (no results file exists)
    with mock.patch("os.path.exists", return_value=False):
        res = run(screen, adata, genes_list, expr_vals, results_file, rank)

    # Ensure that the run function returns the correct format of results
    assert isinstance(res, list)

    # Ensure the screen was called for the batch of experiments
    screen.run_experiment.assert_called()

    # Ensure the pickle file was opened for writing results
    mock_open.assert_called()
    mock_pickle_dump.assert_called()


# Test: run function (resuming from a saved state)


def test_run_resume(mock_open, mock_pickle_dump, mock_comm):
    screen = mock.Mock()  # Mock the PerturbationScreen
    adata = mock.Mock()  # Mock the AnnData object
    genes_list = ["gene1", "gene2", "gene3"]
    expr_vals = [0.5, 1.0, 1.5]
    results_file = "./results.pkl"
    rank = 0

    # Simulate that the results file already exists, and has partial data
    previous_results = [{"experiment": "gene1", "result": 0.5}]
    with mock.patch("os.path.exists", return_value=True):
        with mock.patch("pickle.load", return_value=previous_results):
            res = run(screen, adata, genes_list, expr_vals, results_file, rank)

    # Ensure the results were resumed from the last saved state
    # Mocked experiments should resume from the second experiment
    assert len(res) == 3
    # First experiment should be resumed
    assert res[0]["experiment"] == "gene1"

    # Ensure the pickle file was opened for reading and writing results
    mock_open.assert_called()
    mock_pickle_dump.assert_called()


def test_load_config():
    # Mock configuration data
    config_data = {
        "adata_file": "fib_filtered.h5ad",
        "group_by": "tissue",
        "genes_csv": "gene_list.csv",
        "perturbation_range": [-250, 250],
        "perturbation_steps": 4,
        "perturbation_attempts": 2,
        "screen_type": "DynamoPerturbationScreen",
        "log_file_template": "log-b{rank}.log",
        "results_file_template": "results-b{rank}.pkl",
        "compiled_results_file": "results.pkl",
    }

    # Mock the open function to return our config data as a YAML stream
    with mock.patch("builtins.open", mock.mock_open(read_data="")):
        with mock.patch("yaml.safe_load", return_value=config_data):
            with open("config.yaml") as f:
                config = yaml.safe_load(f)

            # Assert that the config values are correctly loaded
            assert config["adata_file"] == "fib_filtered.h5ad"
            assert config["group_by"] == "tissue"
            assert config["perturbation_range"] == [-250, 250]
            assert config["perturbation_steps"] == 4


def test_import_screen_implementation():
    # Mock the dynamically imported module and class
    mock_module = mock.MagicMock()
    mock_class = mock.MagicMock()

    # Expecting it to import from implementations.dynamo_perturbation_screen
    with mock.patch("importlib.import_module", return_value=mock_module):
        mock_module.DynamoPerturbationScreen = mock_class
        screen_class = import_screen_implementation(
            "DynamoPerturbationScreen", "implementations"
        )

    # Ensure the correct class was retrieved from the module
    assert screen_class == mock_module.DynamoPerturbationScreen
