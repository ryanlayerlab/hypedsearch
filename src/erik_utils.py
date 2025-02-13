import functools
import gzip
import logging
import os
import pickle
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from itertools import chain
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Iterator,
    List,
    Literal,
    Optional,
    Union,
)

import numpy as np
import pandas as pd
from pydantic import BaseModel

from src.erik_constants import MASS, MAX_PEPTIDE_LEN, SEQ
from src.lookups.constants import (
    AMINO_ACID_MASSES,
    DOUBLY_CHARGED_B_BASE,
    DOUBLY_CHARGED_Y_BASE,
    SINGLY_CHARGED_B_BASE,
    SINGLY_CHARGED_Y_BASE,
)


def file_exists(file_name: str) -> bool:
    return os.path.isfile(file_name)


def dataclasses_to_df(dataclasses: List[Any]) -> pd.DataFrame:
    return pd.DataFrame([asdict(item) for item in dataclasses])


def pydantic_models_to_df(models: List[BaseModel]) -> pd.DataFrame:
    return pd.DataFrame([model.model_dump() for model in models])


def make_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def abs_relative_change(obs: Union[float, List[float]], ref: Union[float, List[float]]):
    return np.abs(obs - ref) / np.abs(ref)


def compute_peptide_neutral_mass(
    peptide: str,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    charge: int = 1,
) -> float:
    total_mass = sum([amino_acid_mass_lookup[aa] for aa in peptide])
    return total_mass / charge


def compute_theoretical_mass_over_charge(
    peptide: str,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    charge: int = 1,
) -> float:
    total_mass = sum([amino_acid_mass_lookup[aa] for aa in peptide])
    return total_mass / charge


def setup_logger(
    log_file: str = "app.log",
    mode: Literal["a", "w"] = "a",
    log_level: int = logging.INFO,
):
    """
    Sets up the logger with a specific log file and format.
    """
    # Create a root logger
    logger = logging.getLogger()
    logger.setLevel(log_level)

    # Create a handlers
    file_handler = logging.FileHandler(log_file, mode=mode)
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

    return logger


def log_params(func):
    """
    Decorator to log the parameters of a function when it's called.
    """

    def wrapper(*args, **kwargs):
        # Log function name and its arguments
        logging.info(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


def pickle_object(obj: Any, filename: str):
    """
    Pickles an object with high priority using the latest protocol.

    Args:
        obj: The object to be pickled.
        filename: The name of the file where the pickled object will be stored.
    """
    try:
        with open(filename, "wb") as file:
            # Use the highest protocol for efficiency
            pickle.dump(obj, file, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Object pickled and saved to {filename}")
    except Exception as e:
        print(f"Error while pickling object: {e}")


def unpickle_object(filename: str):
    """
    Unpickles an object with high priority.

    Args:
        filename: The name of the file containing the pickled object.

    Returns:
        The unpickled object.
    """
    try:
        with open(filename, "rb") as file:
            obj = pickle.load(file)
        print(f"Object unpickled from {filename}")
        return obj
    except Exception as e:
        print(f"Error while unpickling object: {e}")
        return None


def pickle_and_compress(obj: Any, file_path: str):
    with gzip.open(file_path, "wb") as file:
        pickle.dump(obj, file, protocol=pickle.HIGHEST_PROTOCOL)


def decompress_and_unpickle(filename: str):
    with gzip.open(filename, "rb") as file:
        return pickle.load(file)


def flatten_list_of_lists(list_of_lists: List[List]):
    return list(chain.from_iterable(list_of_lists))


def run_in_parallel(fcn_of_one_variable: Callable, input_array: List) -> List:
    with ThreadPoolExecutor() as executor:
        result = list(executor.map(fcn_of_one_variable, input_array))

    # Flatten list
    # return flatten_list_of_lists(result)
    return result


def timeit(func):
    """Decorator to time a function and print the execution time."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"Function {func.__name__!r} executed in {elapsed_time:.6f} seconds")
        return result

    return wrapper
