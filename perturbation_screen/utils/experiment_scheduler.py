import logging
from typing import Any, Generator, List, Tuple

import anndata
import numpy as np
from interfaces.perturbation_screen import PerturbationScreen
from utils.logger import log


def range_perturb(
    genes,
    *,
    value_range: Tuple[int, int] = (-1000, 1000),
    steps: int = 17,
    remove_zeros=True,
) -> Tuple[List[List[str]], List[List[int]]]:
    """Performs range perturbations on a list of genes. Returns a list of lists of genes and a list of lists of perturbations."""

    # Check if input is a single gene string
    if isinstance(genes, str):
        genes = [genes]

    # Handle list of lists recursively and zip results
    if any(isinstance(g, list) for g in genes):  # Check if there's any list in genes
        results = [
            range_perturb(sublist, value_range=value_range, steps=steps)
            for sublist in genes
        ]
        combined_genes_list, combined_expr_vals = zip(*results)
        return [item for sublist in combined_genes_list for item in sublist], [
            item for sublist in combined_expr_vals for item in sublist
        ]

    # Handle the regular list of gene strings
    genes_list = [genes] * steps
    expr_vals = [
        [perturb] * len(genes)
        for perturb in np.linspace(*value_range, steps).astype(int)
    ]

    # Remove zeros
    if remove_zeros:
        for i, experiment in enumerate(expr_vals.copy()):
            if all(expr_val for expr_val in experiment) == 0:
                genes_list.pop(i)
                expr_vals.pop(i)

    return genes_list, expr_vals


def batch_split(lst: List, n_splits: int) -> List[List[int]]:
    split_size, remainder = divmod(len(lst), n_splits)
    splits = []
    start = 0

    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)
        splits.append(list(range(start, end)))
        start = end

    return splits


def schedule_experiments(
    genes: List[str], n_partitions: int, **kwargs
) -> Tuple[List[List[str]], List[List[int]], List[List[int]]]:
    # Generate perturbations
    single_gene_experiments = [[g] for g in genes]
    genes_list, expr_vals = range_perturb(single_gene_experiments, **kwargs)

    # Split into batches
    partitions = batch_split(genes_list, n_splits=n_partitions)

    return genes_list, expr_vals, partitions


def batch_experiments(
    screen: PerturbationScreen,
    adata: anndata.AnnData,
    genes_list: List[List[str]],
    expr_vals: List[List[int]],
    *,
    start_idx=0,
    attempts=3,
    logger: logging.Logger | None = None,
    **kwargs,
) -> Generator[dict[str, Any], None, None]:
    # Calculate the total number of experiments as the sum of the current work list and those previously completed
    n_exp = len(genes_list) + start_idx

    genes = np.array(genes_list)
    exprs = np.array(expr_vals)

    for i, (gene, expr) in enumerate(zip(genes, exprs)):
        exp_idx = i + start_idx
        attempt = 0
        success = False
        res = None

        while attempt < attempts and not success:
            try:
                log(
                    f"Running experiment {exp_idx} / {n_exp} | {gene} | {expr} | Attempt {attempt + 1}...",
                    logger=logger,
                )
                res = screen.run_experiment(
                    adata, gene.tolist(), expr.tolist(), **kwargs
                )
                success = True

            except Exception as e:
                log(
                    f"Experiment {exp_idx} failed on attempt {attempt + 1}.",
                    logger=logger,
                    level=logging.ERROR,
                )
                log(str(e), logger=logger, level=logging.ERROR)
                attempt += 1
                # Next attempt with slightly different perturbation
                expr += np.random.randint(-10, 10)

        yield {
            "gene": gene,
            "expr": expr,
            "res": res,
        }
