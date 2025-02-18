import os
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import pandas as pd

GIT_REPO_DIR = Path(__file__).parents[1]

sys.path.append(str(GIT_REPO_DIR))

from src.erik_constants import CHARGE, PRODUCT_ION_TABLE
from src.lookups.protein_product_ion_db import (
    ProteinProductIonDb,
    format_rows_of_product_ion_table,
    get_average_mass_search_time,
    load_existing_protein_product_ion_db,
)


def get_random_sample_of_ions(db: ProteinProductIonDb, sample_size: int):
    ions = db.get_random_sample_of_rows(
        table_name=PRODUCT_ION_TABLE, sample_size=sample_size
    )
    ions = format_rows_of_product_ion_table(rows=ions)
    return ions


@dataclass
class MassSearch:
    mass: float
    time: float
    num_matching_ions: int
    len_search_ion: Optional[int]
    num_matching_ions_by_charge: Counter


def search_over_random_sample_of_ions(
    db: ProteinProductIonDb, sample_size: int, ppm_tolerance: float
) -> List[MassSearch]:
    ions = get_random_sample_of_ions(db=db, sample_size=sample_size)

    search_results = []
    for ion in ions:
        # Perform search
        start_time = time.time()
        matching_ions = db.get_ions_within_mass_tolerance(
            query_mass=ion.mass, ppm_tolerance=ppm_tolerance
        )
        end_time = time.time()

        # Statistics on search
        charges = [ion.charge for ion in matching_ions]
        charge_counts = dict(Counter(charges))
        # charge_counts = {f"{charge}": count for charge, count in Counter(charges).items()}
        result = MassSearch(
            mass=ion.mass,
            time=(end_time - start_time),
            num_matching_ions=len(matching_ions),
            len_search_ion=(ion.exclusive_end - ion.inclusive_start),
            num_matching_ions_by_charge=charge_counts,
        )
        search_results.append(result)
    return search_results


def generate_file_name(db_path: Path, sample_size: int, ppm_tolerance: float):
    return f"{db_path.name}_sampleSize={sample_size}_ppmTol={ppm_tolerance}"


def main(
    db_path: str, sample_size: int, ppm_tolerance: float, output_dir: Optional[str]
):
    db_path = Path(db_path)

    db = load_existing_protein_product_ion_db(db_path=db_path)
    search_results = search_over_random_sample_of_ions(
        db=db, sample_size=sample_size, ppm_tolerance=ppm_tolerance
    )

    df = pd.DataFrame([asdict(result) for result in search_results])

    if output_dir is not None:
        output_dir = Path(output_dir)
        file_name = generate_file_name(
            db_path=db_path, sample_size=sample_size, ppm_tolerance=ppm_tolerance
        )
        output_path = output_dir / f"{file_name}.h5"
        df.to_hdf(output_path, key="df", mode="w", index=False)

    return df


# TESTING
print(GIT_REPO_DIR)
print(os.getcwd())
print(f"Python path: {sys.path}")

# Load database
db_path = "dbs/test.proteinproduction.db"
db = load_existing_protein_product_ion_db(db_path=db_path)
db.print_db_info()

main(db_path=db_path, sample_size=10, ppm_tolerance=10, output_dir="dbs")
