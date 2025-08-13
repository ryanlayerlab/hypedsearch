import gzip
import hashlib
import json
import logging
import os
import pickle
import platform
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from functools import wraps
from itertools import chain
from pathlib import Path
from time import time
from typing import Annotated, Any, Callable, Dict, List, Literal, Optional, Union

import click
import pandas as pd
from pydantic import BeforeValidator
from scipy.stats import percentileofscore

from src.constants import COMET_DIR

MAC_OS = "macos"
LINUX_OS = "linux"

logger = logging.getLogger(__name__)


@dataclass
class CmdLineResult:
    stdout: Optional[str] = None
    stderr: Optional[str] = None
    code: Optional[int] = None


@dataclass
class Position:
    inclusive_start: int
    exclusive_end: int


@dataclass
class Kmer:
    seq: str
    position: Position


def prefixes(seq: str) -> List[str]:
    return [seq[:i] for i in range(1, len(seq) + 1)]


def suffixes(seq: str) -> List[str]:
    return [seq[i:] for i in range(0, len(seq))]


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


def to_path(path: Union[str, Path], check_exists: bool = False) -> Path:
    path = Path(path).absolute()
    if check_exists:
        assert path.exists(), f"The given path {path} does NOT exist!"

    return path


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
        f"{rounding(time_sec / (60 * 60))} hrs"
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


def run_in_parallel(
    fcn_of_one_variable: Callable,
    input_array: List,
    parallel_type: Literal["thread", "process"] = "thread",
) -> List:
    # results = []
    if parallel_type == "thread":
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(fcn_of_one_variable, input_array))
    elif parallel_type == "process":
        with ProcessPoolExecutor() as executor:
            results = list(executor.map(fcn_of_one_variable, input_array))
    else:
        raise ValueError(
            f"Invalid parallel_type: {parallel_type}. Use 'thread' or 'process'."
        )
    return results


def make_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def flatten_list_of_lists(list_of_lists: List[List]):
    return list(chain.from_iterable(list_of_lists))


def mass_difference_in_ppm(mass1: float, mass2: float) -> float:
    return (abs(mass1 - mass2) / mass1) * (10**6)


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


def log_time(level=logging.INFO):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time()
            result = func(*args, **kwargs)
            duration = time() - start
            # logging.debug(f"{func.__name__} took {get_time_in_diff_units(duration)}")
            logging.log(
                level, f"{func.__name__} took {get_time_in_diff_units(duration)}"
            )
            return result

        return wrapper

    return decorator


def run_command_line_cmd(
    cmd: str,
) -> CmdLineResult:
    run_args = shlex.split(cmd)
    result = subprocess.run(run_args, capture_output=True, text=True)
    return CmdLineResult(
        code=result.returncode, stdout=result.stdout, stderr=result.stderr
    )


# function that takes in a list of any one kind dataclass and spits out a dataframe with one row per dataclass list element.
# Optionally take a list of the dataclass's fields/attributes to include in the dataframe.
def dataclass_list_to_df(
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


def get_max_elements(objs: List[Any], attr: str):
    if not objs:
        return []
    max_val = max(getattr(obj, attr) for obj in objs)
    return list(filter(lambda obj: getattr(obj, attr) == max_val, objs))


def get_fcn_of_objects(objs: List[Any], attr: str, fcn: Callable):
    return fcn([getattr(obj, attr) for obj in objs])


def get_arg_fcn_of_objects(objs: List[Any], attr: str, fcn: Callable):
    val = get_fcn_of_objects(objs=objs, attr=attr, fcn=fcn)
    return list(filter(lambda obj: getattr(obj, attr) == val, objs))


def get_os():
    os = platform.platform()
    if "macOS" in os:
        return MAC_OS
    elif "Linux" in os:
        return LINUX_OS
    else:
        raise RuntimeError(f"Unrecognized operating system: {os}")


def get_percentile(value, values):
    """
    E.g., if this function returns 99, it means that 99% of the values in the list
    are less than or equal to value.
    """
    return percentileofscore(values, value, kind="rank")


def get_rank(values, query_val) -> int:
    sorted_unique = sorted(set(values), reverse=True)
    rank_map = {val: rank + 1 for rank, val in enumerate(sorted_unique)}
    return int(rank_map.get(query_val, 0))


def number_greater_than(values: List[float], query_value: float) -> int:
    """
    Returns the number of values in the list that are greater than the query value.
    """
    return sum(1 for val in values if val > query_value)


def check_path_exists(path: Union[Path, str]):
    if isinstance(path, str):
        path = Path(path)
    path = path.absolute()
    if not path.exists():
        raise ValueError(f"Path {path} does not exist")
    return path


ExistingPath = Annotated[Path, BeforeValidator(check_path_exists)]


def get_default_comet_executable_path():
    try:
        os = get_os()
        comet_path = COMET_DIR / f"comet.{os}.exe"
        assert comet_path.exists()
        return comet_path
    except:
        raise RuntimeError(
            "Please provide the Comet executable path. Could NOT find it in the default locations"
        )


def to_json(data: Any, out_path: Union[str, Path]):
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


def load_json(in_path: Union[str, Path]) -> Dict:
    with open(in_path, "r") as f:
        return json.load(f)


def read_new_line_separated_file(path: Union[str, Path]) -> List[str]:
    with open(path, "r") as f:
        lines = [line.strip() for line in f]
        return lines


class PathType(click.ParamType):
    name = "path"

    def convert(self, value, param, ctx):
        try:
            return Path(value).resolve()  # Convert to absolute Path object
        except Exception as e:
            self.fail(f"{value} is not a valid path: {e}", param, ctx)
