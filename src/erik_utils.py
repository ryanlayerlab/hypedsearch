import logging
import os
from dataclasses import dataclass
from typing import Dict, Generator, Iterator, List, Union

import numpy as np
import pandas as pd
from Bio import SeqIO
from pydantic import BaseModel

from src.erik_constants import MASS, MAX_PEPTIDE_LEN, SEQ
from src.lookups.constants import AMINO_ACID_MASSES


def file_exists(file_name: str) -> bool:
    return os.path.isfile(file_name)


def pydantic_models_to_df(models: List[BaseModel]) -> pd.DataFrame:
    return pd.DataFrame([model.model_dump() for model in models])


def make_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def abs_relative_change(obs: Union[float, List[float]], ref: Union[float, List[float]]):
    return np.abs(obs - ref) / np.abs(ref)


def compute_theoretical_mass_over_charge(
    peptide: str,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    charge: int = 1,
):
    total_mass = sum([amino_acid_mass_lookup[aa] for aa in peptide])
    return total_mass / charge


@dataclass
class Kmer:
    seq: str
    mass: float


@dataclass
class Protein:
    # num_in_file: int
    seq: str
    desc: str


def get_proteins_from_fasta(fasta_path: str) -> Iterator[Protein]:
    for protein_num, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
        yield Protein(
            # num_in_file=protein_num,
            seq=str(protein.seq),
            desc=protein.description,
        )


def setup_logger(log_file: str = "app.log", log_level: int = logging.INFO):
    """
    Sets up the logger with a specific log file and format.
    """
    # Create a root logger
    logger = logging.getLogger()
    logger.setLevel(log_level)

    # Create a handlers
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)

    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)

    # Define a common formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to the root logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)


def log_params(func):
    """
    Decorator to log the parameters of a function when it's called.
    """

    def wrapper(*args, **kwargs):
        # Log function name and its arguments
        logging.info(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper
