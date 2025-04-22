import gzip
import hashlib
import logging
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import wraps
from itertools import chain
from pathlib import Path
from time import time
from typing import Any, Callable, List, Literal, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


def remove_gene_name(protein_name: str) -> str:
    """
    Make column for shortened name <database>|<accession number> instead of
    <database>|<accession number>|<gene name>
    """
    return "|".join(protein_name.split("|")[:2])


def make_directory(dir_path: str, overwrite: bool = False) -> None:
    """
    Ensures that a directory exists at the given path.

    Args:
        path (str): The directory path.
        overwrite (bool): If True, deletes the existing directory and creates a new one.
    """
    if os.path.exists(dir_path):
        logger.info(f"{dir_path} already exists")
        if overwrite:
            logger.info(f"Overwriting {dir_path}")
            shutil.rmtree(dir_path)  # Remove the existing directory and its contents
            os.makedirs(dir_path)  # Create a new empty directory
    else:
        logger.info(f"{dir_path} does not exist. Creating it...")
        os.makedirs(dir_path)  # Create the directory if it does not exist


def to_path(path: Union[str, Path]):
    return Path(path)


@dataclass
class Position:
    inclusive_start: int
    exclusive_end: int


@dataclass
class Kmer:
    seq: str
    position: Position

    # def as_product_ion(self, charge: int, ion_type: IonTypes, )


def generate_aa_kmers(aa_seq: str, max_k: int, min_k: Optional[int] = 1) -> List[Kmer]:
    disallowed_amino_acid_symbols = ["B", "X", "U", "Z", "O", "J"]
    kmers = []
    for kmer_len in range(min_k, max_k + 1):  # Iterate over all k from 1 to max_k
        for inclusive_start in range(
            len(aa_seq) - kmer_len + 1
        ):  # Slide over the string
            exclusive_end = inclusive_start + kmer_len
            kmer = aa_seq[inclusive_start:exclusive_end]
            if any(char in disallowed_amino_acid_symbols for char in kmer):
                continue
            kmers.append(
                Kmer(
                    seq=kmer,
                    position=Position(
                        inclusive_start=inclusive_start,
                        exclusive_end=exclusive_end,
                    ),
                )
            )

    return kmers


def get_positions_of_subseq_in_seq(subseq: str, seq: str) -> List[Position]:
    matches = re.finditer(subseq, seq)
    return [
        Position(inclusive_start=match.start(), exclusive_end=match.end())
        for match in matches
    ]


def setup_logger(
    log_file: str = "app.log",
    mode: Literal["a", "w"] = "a",
    log_level: int = logging.INFO,
) -> logging.Logger:
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
    console_handler = logging.StreamHandler(stream=sys.stdout)
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


def get_time_in_diff_units(time_sec: float, decimal_places: int = 2) -> str:
    rounding = lambda x: round(x, decimal_places)
    return (
        f"{rounding(time_sec)} secs = "
        f"{rounding(time_sec / 60)} mins = "
        f"{rounding(time_sec / (60*60))} hrs"
    )


def log_params(func):
    """
    Decorator to log the parameters of a function when it's called.
    """

    def wrapper(*args, **kwargs):
        # Log function name and its arguments
        logging.info(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        return func(*args, **kwargs)

    return wrapper


def relative_ppm_tolerance_in_daltons(ppm: float, ref_mass: float) -> float:
    return (ppm * ref_mass) / (10**6)


def run_in_parallel(fcn_of_one_variable: Callable, input_array: List) -> List:
    # results = []
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(fcn_of_one_variable, input_array))

        # futures = {
        #     executor.submit(fcn_of_one_variable, item): item for item in input_array
        # }

        # # Use tqdm to show progress
        # for future in tqdm(
        #     as_completed(futures), total=len(input_array), desc="Processing"
        # ):
        #     results.append(future.result())
    # Flatten list
    # return flatten_list_of_lists(result)
    return results


def make_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def flatten_list_of_lists(list_of_lists: List[List]):
    return list(chain.from_iterable(list_of_lists))


def mass_difference_in_ppm(ref_mass: float, query_mass: float) -> float:
    return (abs(ref_mass - query_mass) / ref_mass) * (10**6)


def pickle_and_compress(obj: Any, file_path: str):
    with gzip.open(file_path, "wb") as file:
        pickle.dump(obj, file, protocol=pickle.HIGHEST_PROTOCOL)


def decompress_and_depickle(file_path: str):
    with gzip.open(file_path, "rb") as file:
        return pickle.load(file)


def file_hash(filepath: Path, algorithm="sha256") -> str:
    """Compute the hash of a file using the specified algorithm."""
    hash_func = hashlib.new(algorithm)
    with filepath.open("rb") as f:
        for chunk in iter(
            lambda: f.read(4096), b""
        ):  # Read in chunks to handle large files efficiently
            hash_func.update(chunk)
    return hash_func.hexdigest()


def log_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        duration = time() - start
        logging.info(f"{func.__name__} took {get_time_in_diff_units(duration)}")
        return result

    return wrapper


@dataclass
class CmdResult:
    stdout: Optional[str] = None
    stderr: Optional[str] = None
    code: Optional[int] = None


def run_command_line_cmd(cmd: str):
    run_args = shlex.split(cmd)
    result = subprocess.run(run_args, capture_output=True, text=True)
    return CmdResult(code=result.returncode, stdout=result.stdout, stderr=result.stderr)


# function that takes in a list of any one kind dataclass and spits out a dataframe with one row per dataclass list element.
# Optionally take a list of the dataclass's fields/attributes to include in the dataframe.
def dataclass_list_to_dataframe(
    dataclass_list: List[Any], fields: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Function that takes in a list of any one kind dataclass and spits out a dataframe
    with one row per dataclass list element. Optionally take a list of the dataclass's
    fields/attributes to include in the dataframe.
    """
    if fields is None:
        fields = [
            field.name for field in dataclass_list[0].__dataclass_fields__.values()
        ]
    data = {field: [] for field in fields}
    for item in dataclass_list:
        for field in fields:
            data[field].append(getattr(item, field))
    return pd.DataFrame(data)


def lowercase_and_underscore(in_str: str):
    """
    Change a string to all lower case and replace spaces with underscores.
    """
    return in_str.lower().replace(" ", "_")


def is_pickleable(obj: Any) -> bool:
    """
    Check if an object is pickleable.
    """
    try:
        pickle.dumps(obj)
        return True
    except Exception as e:
        print(f"Not pickleable: {e}")
        return False
