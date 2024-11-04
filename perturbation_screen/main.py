import logging
import os
import pickle
import warnings
from typing import Any, List

import anndata
import numpy as np
from interfaces.perturbation_screen import (
    PerturbationScreen,
    filter_genes,
    import_screen_implementation,
)
from mpi4py import MPI
from utils.config import load_config
from utils.data_promise import DataPromise
from utils.experiment_scheduler import batch_experiments, schedule_experiments
from utils.logger import log, setup_logger

warnings.filterwarnings("ignore")


def save_results(res, path):
    with open(path, "wb") as f:
        pickle.dump(res, f)


def main():
    # Load the config
    config = load_config("config.yaml")

    # Initialize the MPI environment
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    os.makedirs(config.working_directory, exist_ok=True)
    os.chdir(config.working_directory)

    # Parameters
    log_file = config.log_file_template.format(rank=rank, size=size)
    results_file = config.results_file_template.format(rank=rank, size=size)
    group_label = config.group_by

    # Instantiate the logger
    logger: logging.Logger = setup_logger("PerturbationExperiment", log_file)

    # Instantiate the perturbation screen
    PerturbationScreenImplementation = import_screen_implementation(config.screen_type)
    screen: PerturbationScreen = PerturbationScreenImplementation(
        group_label=group_label, logger=logger
    )

    # Create the DataPromise
    data_promise = DataPromise(config, logger=logger)

    # Load the data
    adata = data_promise.resolve(screen, comm)

    if rank == 0:
        # Retrieve group labels
        group_labels = screen.get_group_labels(adata)
        # Load genes
        with open(config.genes_csv, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        log(f"Loaded {len(genes)} genes from {config.genes_csv}.", logger=logger)
        # Filter genes
        included_genes, excluded_genes = filter_genes(
            genes, screen.get_valid_genes(adata, gene_list=genes)
        )
        log(
            f"Excluded {len(excluded_genes)} genes from {len(genes)} total.",
            logger=logger,
        )
        if excluded_genes:
            log(f"Excluded genes: {excluded_genes}", logger=logger)
        # Schedule perturbations
        genes_list, expr_vals, partitions = schedule_experiments(
            genes=included_genes,
            n_partitions=size,
            steps=config.perturbation_steps,
            value_range=config.perturbation_range,
        )
    else:
        group_labels = []
        included_genes, excluded_genes = [], []
        genes_list, expr_vals, partitions = [], [], []

    # Broadcast data to all processes
    (
        group_labels,
        included_genes,
        excluded_genes,
        genes_list,
        expr_vals,
        partitions,
    ) = comm.bcast(
        (
            group_labels,
            included_genes,
            excluded_genes,
            genes_list,
            expr_vals,
            partitions,
        ),
        root=0,
    )
    log(f"Broadcasted/received ancillary data successfully.", logger=logger)

    # Run the screen
    log(f"Running screen on partition {rank + 1}/{size}.", logger=logger)
    partition = partitions[rank]
    res = run(
        screen,
        adata,
        genes_list=np.array(genes_list)[partition],
        expr_vals=np.array(expr_vals)[partition],
        results_file=results_file,
        rank=rank,
        logger=logger,
    )

    log(f"Parition {rank + 1}/{size} screen complete.", logger=logger)

    results = comm.gather(res, root=0)

    if rank == 0 and results:
        export_dict = {
            "included_genes": included_genes,
            "excluded_genes": excluded_genes,
            "group_label": group_label,
            "group_labels": group_labels,
            "results": sum(results, []),
            "perturbation_range": config.perturbation_range,
            "perturbation_steps": config.perturbation_steps,
            "perturbation_attempts": config.perturbation_attempts,
            "screen_type": config.screen_type,
        }
        save_results(export_dict, config.compiled_results_file)
        log(
            f"Gathered results saved to `{config.compiled_results_file}`.",
            logger=logger,
        )


def run(
    screen: PerturbationScreen,
    adata: anndata.AnnData,
    genes_list,
    expr_vals,
    results_file,
    rank,
    logger=None,
    **kwargs,
):
    start_idx = 0
    res: List[dict[str, Any]] = []

    # Check to see if we are starting mid-way through a run
    if os.path.exists(results_file):
        with open(results_file, "rb") as f:
            res = pickle.load(f)
        start_idx = len(res)
        log(f"Resuming screen (b{rank}) from experiment {start_idx}...", logger=logger)

    else:
        log(f"Starting new screen (b{rank}).", logger=logger)

    # loop through experiments, saving as we go
    for exp in batch_experiments(
        screen,
        adata,
        genes_list[start_idx:],
        expr_vals[start_idx:],
        start_idx=start_idx,
        logger=logger,
        **kwargs,
    ):
        res.append(exp)
        save_results(res, path=results_file)

    return res


if __name__ == "__main__":
    main()
