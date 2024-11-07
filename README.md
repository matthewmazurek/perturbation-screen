# Perturbation Screen Pipeline

This repository contains a Python implementation of a parallelized perturbation screening pipeline designed for single-cell RNA-seq data analysis. The pipeline leverages MPI for distributed computing, allowing efficient processing of large datasets.

## Features

- **Parallel Processing with MPI**: Utilize multiple processors to accelerate computations using `mpi4py`.
- **Customizable Perturbation Screens**: Implement your own perturbation screens by extending the `PerturbationScreen` interface.
- **Configurable Settings**: Adjust parameters easily via a YAML configuration file.
- **Robust Logging and Checkpointing**: Keep track of progress and resume computations with logging and result persistence.

## Requirements

- **Python 3.10**
- [NumPy](https://numpy.org/)
- [AnnData](https://anndata.readthedocs.io/en/latest/)
- [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
- [PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation) (for YAML configuration)
- [Dynamo](https://dynamo-release.readthedocs.io/en/latest/) (for Dynamo implementation)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/matthewmazurek/perturbation-screen.git
cd perturbation-screen
```

### 2. Install or Activate MPI

MPI is required for parallel processing. Install MPI or activate it if it's already available.

- **On Cloud Computing Systems or HPC Clusters with MPI pre-installed**:

  Activate the OpenMPI module using your system's module environment. For example:

  ```bash
  module load openmpi
  ```

  You can list available OpenMPI modules with:

  ```bash
  module avail openmpi
  ```

- **Ubuntu**:

  ```bash
  sudo apt-get install libopenmpi-dev openmpi-bin
  ```

- **macOS** (using Homebrew):

  ```bash
  brew install open-mpi
  ```

---

### 3. Set up the Conda Environment

It's recommended to use Conda for managing dependencies and ensuring compatibility with Python 3.10.

1. **Create a new Conda environment**:

   ```bash
   conda create -n pspipe python=3.10
   ```

2. **Activate the environment**:

   ```bash
   conda activate pspipe
   ```

3. **Install the required packages**:

   ```bash
   pip install -r requirements.txt
   ```

## Configuration

The pipeline is configured via a `config.yaml` file located in the project's root directory. Below is an example of the configuration options:

```yaml
working_directory: './results'
log_file_template: 'log-b{rank}.log'
results_file_template: 'results-b{rank}.pkl'
compiled_results_file: 'compiled_results.pkl'
group_by: 'cell_type'
screen_type: 'DynamoPerturbationScreen'
genes_csv: 'genes_list.csv'
perturbation_range: [-1000, 1000]
perturbation_steps: 11
perturbation_attempts: 3
```

- **working_directory**: Directory where results and logs will be stored.

- **log_file_template**: Template for log file naming. You can use the following template variables:

  - `{rank}`: The rank of the MPI process.
  - `{size}`: The total number of MPI processes.

- **results_file_template**: Template for individual results file naming. Supports the same template variables as `log_file_template`.

- **compiled_results_file**: Name of the file where compiled results from all processes are saved.

- **group_by**: The key in the AnnData object to group cells by.

- **screen_type**: The type of perturbation screen to use. Available options:

  - `'DynamoPerturbationScreen'`: Uses Dynamo for perturbation analysis.
  - `'DummyPerturbationScreen'`: A dummy implementation for testing purposes.

- **genes_csv**: Path to a CSV file containing a list of genes to perturb.

- **perturbation_range**: Range of expression values to test for perturbations.

- **perturbation_steps**: Number of equidistant points to sample within the perturbation range.

- **perturbation_attempts**: Number of attempts for each perturbation before passing.

## Usage

### 1. Prepare your data

- Ensure you have an `AnnData` object containing your single-cell RNA-seq data.
- Create a CSV file (`genes_list.csv`) with the list of genes you want to perturb.

### 2. Configure the pipeline

- Edit the `config.yaml` file to match your data and desired parameters.

### 3. Run the pipeline

Execute the main script using `srun` to utilize MPI for parallel processing. The `--mpi=pmi2` flag ensures compatibility with certain MPI implementations.

```bash
srun --mpi=pmi2 -n [number_of_processes] python perturbation_screen/main.py
```

**Example**:

```bash
srun --mpi=pmi2 -n 4 python perturbation_screen/main.py
```

### 4. Monitor the progress

- Logs will be saved in the `working_directory` as specified in the configuration.
- Individual results are saved incrementally, allowing you to resume computation if interrupted.

### 5. Results

- Upon completion, compiled results will be saved to the file specified by `compiled_results_file` in the configuration.
- Use appropriate tools to analyze and visualize the results.

## Customization

### Implementing Custom Perturbation Screens

To create a custom perturbation screen, you need to extend the `PerturbationScreen` interface. See `DynamoPerturbationScreen`

#### Steps to Implement a Custom Screen

1. **Create a new class** in the `interfaces/perturbation_screen` module that inherits from `PerturbationScreen`.

2. **Implement the required methods**:
   - `load_data(self, adata_file: str) -> anndata.AnnData`
   - `preprocess(self, adata: anndata.AnnData) -> None`
   - `get_group_labels(self, adata: anndata.AnnData) -> List[str]`
   - `get_valid_genes(self, adata: anndata.AnnData, gene_list: List[str] | None = None) -> List[str]`
   - `run_experiment(self, adata: anndata.AnnData, genes: List[str], expr_vals: List[int]) -> Any`

3. **Update the `screen_type`** in `config.yaml` to your new class name.

#### Example: DynamoPerturbationScreen

An example implementation using Dynamo for perturbation analysis is provided as `DynamoPerturbationScreen`. Below is a simplified version:

```python
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
```

### Modifying the Pipeline

Feel free to modify the main script (`main.py`) or any of the utility modules to suit your specific needs.

## Contributing

Contributions are welcome! Please follow these steps:

1. **Fork the repository**.
2. **Create a new branch** for your feature or bugfix.
3. **Make your changes** and commit them with clear messages.
4. **Submit a pull request** detailing your changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE.txt) file for details.

## Contact

For questions or support, please open an issue in the repository.