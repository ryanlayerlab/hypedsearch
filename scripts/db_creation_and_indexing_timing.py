import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

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
from src.peptides_and_ions import Peptide
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


logger = setup_logger()

peptides = Peptide.from_fasta(fasta_path=MOUSE_PROTEOME)

testing = True
max_num_mins = 60
num_peptides_step = 5
num_repeats = 5
charges, ion_types = [1, 2, 3], [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]


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
)

logger.info(f"Running with params = {params}")
data = []
prev_time_delta = 0
idx = 1
while prev_time_delta < 60 * max_num_mins:
    num_peptides = idx * num_peptides_step
    logger.info(f"|peptides| = {num_peptides}; |time delta| = {prev_time_delta}")
    for _ in range(num_repeats):
        # Sample peptides
        peptide_sample = random.sample(peptides, num_peptides)

        # Create and index database
        t0 = time()
        db = create_and_populate_protein_and_product_ion_database(
            charges=charges,
            ion_types=ion_types,
            peptides=peptide_sample,
        )
        db.add_index(
            table_name=PRODUCT_ION_TABLE,
            index_name="mass_idx",
            colms_to_index=[NEUTRAL_MASS],
        )
        prev_time_delta = time() - t0

        # Validation
        assert len(db.get_proteins()) == num_peptides

        # Updates
        data.append([num_peptides, prev_time_delta, len(db.get_product_ions())])
    idx += 1

data = pd.DataFrame(
    data, columns=["num_peptides", "create_and_index_time_delta", "num_product_ions"]
)

file_name = f"db_creation_timing_charge={charges}_maxMins={max_num_mins}.pkl"
file_path = RESULTS_DIR / file_name
if testing:
    file_path = RESULTS_DIR / ("TESTING_" + file_name)
pickle_and_compress(
    obj=data,
    file_path=file_path,
)
