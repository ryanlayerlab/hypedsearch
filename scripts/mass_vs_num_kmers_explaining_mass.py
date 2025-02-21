import time
from dataclasses import dataclass
from typing import Any, List, Tuple

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes

from src.constants import MASS, PLOTS_DIR, PRODUCT_ION_TABLE, SEQ
from src.erik import (
    get_query_to_select_rows_by_mass,
    query_database,
    relative_ppm_tolerance_in_daltons,
)
from src.erik_utils import run_in_parallel
from src.lookups.data_classes import Kmer, KmerWithMass
from src.plot_utils import fig_setup, finalize, set_title_axes_labels

KMERS_DB_PATH = "./dbs/Uniprot_mouse.fasta.uniq_kmers.max_k=50.db"


def sample_from_kmer_db(
    num_rows: int, db_path: str = KMERS_DB_PATH, table_name: str = PRODUCT_ION_TABLE
) -> KmerWithMass:
    query = f"""
    SELECT {SEQ}, {MASS}
    FROM {table_name}
    ORDER BY RANDOM()
    LIMIT {num_rows};
    """
    t0 = time.time()
    rows = query_database(query=query, db_path=db_path)
    t1 = time.time()
    print(f"get_random_rows with num_rows={num_rows} took {round(t1-t0, 2)} seconds")
    kmers = [KmerWithMass(kmer=Kmer(seq=row[0]), mass=row[1]) for row in rows]
    return kmers


def get_num_of_explanatory_kmers(mass: float, ppm_tol: int, db_path: str):
    da_tol = relative_ppm_tolerance_in_daltons(ppm=ppm_tol, ref_mass=mass)
    query = get_query_to_select_rows_by_mass(
        mass=mass,
        mass_tol=da_tol,
        # query_type="both"
    )
    result = query_database(query=query, db_path=db_path)
    num_explanatory_kmers = result[0][0]
    # print(f"mass={mass}\nda_tol={da_tol}\nquery:\n{query}")
    return num_explanatory_kmers


@dataclass
class ExplanatoryKmers:
    query_mass: float
    num_explanatory_kmers: int


def get_explanatory_kmers_for_many_masses(
    masses: List[float], ppm_tol: int, db_path: str
) -> List[ExplanatoryKmers]:
    return [
        ExplanatoryKmers(
            query_mass=mass,
            num_explanatory_kmers=get_num_of_explanatory_kmers(
                mass=mass, ppm_tol=ppm_tol, db_path=db_path
            ),
        )
        for mass in masses
    ]


def sample_masses_and_get_num_explanatory_kmers(
    num_samples: int, kmer_db_path: str, ppm_tol: float
) -> List[ExplanatoryKmers]:
    kmers = sample_from_kmer_db(num_rows=num_samples, db_path=kmer_db_path)
    masses = [kmer.mass for kmer in kmers]
    return get_explanatory_kmers_for_many_masses(
        masses=masses, ppm_tol=ppm_tol, db_path=kmer_db_path
    )


def sample_in_parallel(
    num_samples: int, num_folds: int, ppm_tol: float, kmer_db_path: str
) -> List[ExplanatoryKmers]:
    fcn = lambda num_samples: sample_masses_and_get_num_explanatory_kmers(
        num_samples=num_samples,
        kmer_db_path=kmer_db_path,
        ppm_tol=ppm_tol,
    )
    input_array = [num_samples] * num_folds
    explanatory_kmers = run_in_parallel(
        fcn_of_one_variable=fcn, input_array=input_array
    )
    return explanatory_kmers


def plot_explanatory_kmers(
    # ax: Axes,
    explanatory_kmers: List[ExplanatoryKmers],
    save_path: str,
    ppm_tol: int,
):
    fig, axs = fig_setup()
    fig_setup
    data = np.array(
        list(
            [kmer.query_mass, kmer.num_explanatory_kmers] for kmer in explanatory_kmers
        )
    )
    axs[0].plot(data[:, 0], data[:, 1], ".", ms=0.5)
    set_title_axes_labels(
        axs[0],
        title=f"PPM tolerance = {ppm_tol}",
        xlabel="Mass",
        ylabel="Number explanatory kmers",
    )
    finalize(axs)
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
