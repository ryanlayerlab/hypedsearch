import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

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
from src.utils import pickle_and_compress, setup_logger


@dataclass
class ExperimentParams:
    testing: bool
    max_num_mins: int
    num_peptides_step: int
    num_repeats: int
    charges: List[int]
    ion_types: List[IonTypes]
    kmer_size: int
    num_mass_queries: int


@dataclass
class ExperimentResult:
    params: ExperimentParams
    df: pd.DataFrame


def get_experiment_results_file_path(experiment_params: ExperimentParams):
    # Get operating system
    if "macOS" in platform.platform():
        operating_sys = "MAC"
    # Running on Fiji
    else:
        operating_sys = "FIJI"
    file_name = (
        f"db_timing_os={operating_sys}_" +
        f"charges={experiment_params.charges}_" + 
        f"maxMins={experiment_params.max_num_mins}_" + 
        f"peptideStep={experiment_params.num_peptides_step}_" +
        f"numRepeats={experiment_params.num_repeats}_" + 
        ".pkl"
    )
    file_path = RESULTS_DIR / file_name
    if experiment_params.testing:
        file_path = RESULTS_DIR / ("TESTING_" + file_name)
    return file_path


def db_creation_indexing_querying_timing(experiment_params: ExperimentParams):
    logger.info(f"Running with params = {experiment_params}")

    # Load peptides
    peptides = Peptide.from_fasta(fasta_path=MOUSE_PROTEOME)

    # Experiment
    data = []
    creation_time = 0
    idx = 1
    while creation_time < 60 * experiment_params.max_num_mins:
        num_peptides = idx * experiment_params.num_peptides_step
        logger.info(f"|peptides| = {num_peptides}; |time delta| = {creation_time}")
        for _ in range(num_repeats):
            # Sample peptides
            peptide_sample = random.sample(peptides, num_peptides)

            # Create database
            t0 = time()
            db = create_and_populate_protein_and_product_ion_database(
                charges=experiment_params.charges,
                ion_types=experiment_params.ion_types,
                peptides=peptide_sample,
            )
            creation_time = time() - t0
            # Validation
            assert len(db.get_proteins()) == num_peptides

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
                k=experiment_params.kmer_size,
                sample_size=experiment_params.num_mass_queries,
                peptides=peptide_sample,
            )
            for kmer in query_kmers:
                charge = random.choice(charges)
                ion_type = random.choice(ion_types)
                ion = ProductIon(seq=kmer, charge=charge, ion_type=ion_type)
                t0 = time()
                db.get_ions_within_mass_tolerance(
                    query_mass=ion.neutral_mass, ppm_tolerance=10
                )
                query_times.append(time() - t0)
            avg_query_time = np.mean(query_times)

            # Updates
            data.append(
                [
                    num_peptides,
                    creation_time,
                    index_time,
                    avg_query_time,
                    len(db.get_product_ions()),
                ]
            )

        idx += 1

    data = pd.DataFrame(
        data,
        columns=[
            "num_peptides",
            "create_time",
            "index_time",
            "avg_query_time",
            "num_product_ions",
        ],
    )

    result = ExperimentResult(params=experiment_params, df=data)
    file_path = get_experiment_results_file_path(experiment_params=experiment_params)
    logger.info(f"Pickling data to path {file_path}")
    pickle_and_compress(
        obj=result,
        file_path=file_path,
    )


logger = setup_logger()

max_num_mins = 60
num_peptides_step = 75
num_repeats = 3
charges, ion_types = [1, 2, 3], [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
num_mass_queries = 5
kmer_size = 15


testing = False
if testing:
    max_num_mins = 0.1
    num_peptides_step = 10
    num_repeats = 1

params = ExperimentParams(
    testing=testing,
    max_num_mins=max_num_mins,
    num_peptides_step=num_peptides_step,
    num_repeats=num_repeats,
    charges=charges,
    ion_types=ion_types,
    kmer_size=kmer_size,
    num_mass_queries=num_mass_queries,
)

db_creation_indexing_querying_timing(experiment_params=params)
