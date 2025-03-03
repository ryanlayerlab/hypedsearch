import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

import click
import numpy as np

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

import random
from time import time

import pandas as pd

from src.constants import (
    MOUSE_PROTEOME,
    NEUTRAL_MASS,
    PRODUCT_ION_TABLE,
    RESULTS_DIR,
    IonTypes,
)
from src.peptides_and_ions import Peptide, ProductIon, random_sample_of_unique_kmers
from src.protein_product_ion_database import (
    create_and_populate_protein_and_product_ion_database,
)
from src.utils import make_directory, pickle_and_compress, setup_logger

logger = setup_logger()


@dataclass
class ExperimentParams:
    testing: bool
    # max_num_mins: int
    # num_peptides_step: int
    num_peptides_to_add: int
    num_repeats: int
    charges: List[int]
    ion_types: List[IonTypes]
    query_kmer_size: int
    num_mass_queries: int
    fasta_path: str
    ppm_tolerance: float


@dataclass
class ExperimentResult:
    params: ExperimentParams
    df: pd.DataFrame


def db_creation_indexing_querying_timing(
    experiment_params: ExperimentParams, output_dir: Path
):
    logger.info(f"Running with params = {experiment_params}")

    # Load peptides
    peptides = Peptide.from_fasta(fasta_path=experiment_params.fasta_path)

    # Experiment
    data = []
    # creation_time = 0
    # idx = 1
    # while creation_time < 60 * experiment_params.max_num_mins:
    # num_peptides = idx * experiment_params.num_peptides_step
    # logger.info(f"|peptides| = {num_peptides}; |time delta| = {creation_time}")
    for _ in range(experiment_params.num_repeats):
        # Sample peptides
        peptide_sample = random.sample(peptides, experiment_params.num_peptides_to_add)

        # Create database
        t0 = time()
        db = create_and_populate_protein_and_product_ion_database(
            charges=experiment_params.charges,
            ion_types=experiment_params.ion_types,
            peptides=peptide_sample,
        )
        creation_time = time() - t0
        # Validation
        assert len(db.get_proteins()) == experiment_params.num_peptides_to_add

        # Index database
        t0 = time()
        db.add_index(
            table_name=PRODUCT_ION_TABLE,
            index_name="mass_idx",
            colms_to_index=[NEUTRAL_MASS],
        )
        index_time = time() - t0

        # Query database
        query_times = []
        query_kmers = random_sample_of_unique_kmers(
            k=experiment_params.query_kmer_size,
            sample_size=experiment_params.num_mass_queries,
            peptides=peptide_sample,
        )
        for kmer in query_kmers:
            charge = random.choice(experiment_params.charges)
            ion_type = random.choice(experiment_params.ion_types)
            ion = ProductIon(seq=kmer, charge=charge, ion_type=ion_type)
            t0 = time()
            db.get_ions_within_mass_tolerance(
                query_mass=ion.neutral_mass,
                ppm_tolerance=experiment_params.ppm_tolerance,
            )
            query_times.append(time() - t0)
        avg_query_time = np.mean(query_times)
        sd_query_time = np.std(query_times)

        # Updates
        data.append(
            [
                experiment_params.num_peptides_to_add,
                creation_time,
                index_time,
                avg_query_time,
                sd_query_time,
                len(db.get_product_ions()),
            ]
        )

    data = pd.DataFrame(
        data,
        columns=[
            "num_peptides",
            "create_time",
            "index_time",
            "avg_query_time",
            "sd_query_time",
            "num_product_ions",
        ],
    )

    result = ExperimentResult(params=experiment_params, df=data)
    output_path = output_dir / f"{experiment_params.num_peptides_to_add}.pkl"
    logger.info(f"Pickling data to path {output_path}")
    pickle_and_compress(
        obj=result,
        file_path=output_path,
    )


def get_experiment_results_directory(experiment_params: ExperimentParams) -> Path:
    # Get operating system
    if "macOS" in platform.platform():
        operating_sys = "MAC"
    # Running on Fiji
    else:
        operating_sys = "FIJI"

    dir_name = (
        f"db_timing_os={operating_sys}_"
        + f"charges={experiment_params.charges}_"
        + f"queryKmerSize={experiment_params.query_kmer_size}_"
        + f"numRepeats={experiment_params.num_repeats}"
    )
    # file_name = (
    #     f"db_timing_os={operating_sys}_"
    #     + f"charges={experiment_params.charges}_"
    #     + f"maxMins={experiment_params.max_num_mins}_"
    #     + f"numRepeats={experiment_params.num_repeats}_"
    #     + ".pkl"
    # )

    if experiment_params.testing:
        return RESULTS_DIR / ("TESTING_" + dir_name)
    else:
        return RESULTS_DIR / dir_name


@click.command(
    help="""
    Create protein-product ion database(s)\n
    """,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--num_peptides_to_add",
    "-n",
    type=int,
    required=True,
    # help="Path to FASTA file"
)
@click.option(
    "--testing",
    "-t",
    required=False,
    is_flag=True,
    show_default=True,
    default=False,
    help="Whether to run in testing mode.",
)
@click.option(
    "--num_repeats",
    "-r",
    type=int,
    required=False,
    default=1,
    show_default=True,
    # help="Path to FASTA file"
)
# @click.option(
#     "--charges",
#     "-c",
#     type=int,
#     multiple=True,
#     required=False,
#     help=(
#         "Integer charges to consider. Need to be provided with a '-c' before each charge, "
#         "e.g., -c 1 -c 2 -c 3."
#     ),
# )
def main(
    num_peptides_to_add: int,
    testing: bool = False,
    # max_num_mins: int = 60,
    num_repeats: int = 1,
    charges: List[int] = [1, 2, 3],
    ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
    query_kmer_size: int = 15,
    fasta_path: str = MOUSE_PROTEOME,
    num_mass_queries: int = 100,
    ppm_tolerance: float = 10,
):
    # Set experiment parameters
    params = ExperimentParams(
        testing=testing,
        num_peptides_to_add=num_peptides_to_add,
        # max_num_mins=max_num_mins,
        num_repeats=num_repeats,
        charges=charges,
        ion_types=ion_types,
        query_kmer_size=query_kmer_size,
        fasta_path=fasta_path,
        num_mass_queries=num_mass_queries,
        ppm_tolerance=ppm_tolerance,
    )

    # Get output directory
    output_dir = get_experiment_results_directory(experiment_params=params)
    make_directory(dir_path=output_dir)

    db_creation_indexing_querying_timing(
        experiment_params=params, output_dir=output_dir
    )
    pass


if __name__ == "__main__":
    main()
