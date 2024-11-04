import os

import anndata
from interfaces.perturbation_screen import PerturbationScreen
from mpi4py import MPI
from utils.config import ConfigModel
from utils.logger import log


class DataPromise:
    def __init__(self, config: ConfigModel, load_processed=False, logger=None):
        """
        Initialize a DataPromise object.

        Parameters:
        - config (ConfigModel): User-specified configuration object.
        - load_processed (bool): Flag indicating whether to load processed data if available.
        """
        self.adata = None
        self.adata_file = config.adata_file
        self.adata_file_processed = config.adata_file_processed
        self.load_processed = config.load_processed
        self.logger = logger

    def resolve(
        self, perturbation_screen: PerturbationScreen, comm: MPI.Comm
    ) -> anndata.AnnData:
        """
        Load the data, performing preprocessing if necessary.

        Parameters:
        - perturbation_screen: An instance of PerturbationScreen or its subclass.
        - comm: MPI communicator.

        Returns:
        - adata: The loaded and preprocessed data.
        """
        log(f"Attempting to resolve data promise...", logger=self.logger)
        if self.adata is not None:
            log(f"Data is already in memory.", logger=self.logger)
            return self.adata

        rank = comm.Get_rank()

        # Try to load the processed data if the flag is set and file exists
        if self.load_processed and os.path.exists(self.adata_file_processed):
            self.adata = perturbation_screen.load_data(self.adata_file_processed)
            log(f"Loaded processed data successfully.", logger=self.logger)
            return self.adata

        # Data requires preprocessing
        if rank == 0:
            # Rank 0 loads and preprocesses the data
            self.adata = perturbation_screen.load_data(self.adata_file)
            perturbation_screen.preprocess(self.adata)
            log(f"Preprocessed data successfully.", logger=self.logger)
            if self.load_processed:
                # Save the preprocessed data for future use
                perturbation_screen.save_data(self.adata, self.adata_file_processed)
        else:
            # Other ranks set adata to None
            self.adata = None

        # Synchronize all processes
        comm.Barrier()

        # After synchronization, other ranks load the processed data if saved
        if self.load_processed and os.path.exists(self.adata_file_processed):
            if rank != 0:
                self.adata = perturbation_screen.load_data(self.adata_file_processed)
                log(f"Loaded processed data successfully.", logger=self.logger)
            return self.adata  # type: ignore

        # Rank 0 broadcasts the data if not saved to disk
        else:
            self.adata = comm.bcast(self.adata, root=0)
            log(f"Broadcasted/received data successfully.", logger=self.logger)

        return self.adata
