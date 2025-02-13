import logging
import sys
import time
from pathlib import Path
from typing import List

import click

from src.erik_constants import MAX_PEPTIDE_LEN
from src.erik_utils import log_params, setup_logger
from src.lookups.protein_product_ion_db import (
    ProteinProductIonDb,
    create_protein_product_ion_db,
)

logger = setup_logger()


def database_protein_product_ion_db_filename(
    fasta_path: Path, max_k: int, charges_to_consider: List[int]
):
    return f"{fasta_path.name}_max_k={max_k}_charges={charges_to_consider}.db"


@click.command(
    help="""
    This script creates a database with a table with colums:\n
    <kmer seq>, <kmer mass>, <charge>\n
    """,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--fasta_path", "-f", type=str, required=True, help="Path to FASTA file")
@click.option(
    "--db_dir",
    "-d",
    type=str,
    required=True,
    help="Folder where you'd like the database file to be created.",
)
@click.option(
    "--max_kmer_len",
    "-k",
    help=f"Maximum kmer length to add to database.",
    type=int,
    required=True,
)
@click.option(
    "--charges",
    "-c",
    type=int,
    multiple=True,
    required=True,
    help=(
        "Integer charges to consider. Need to be provided with a '-c' before each charge, "
        "e.g., -c 1 -c 2 -c 3."
    ),
)
@log_params
def create_protein_and_product_ion_database(
    fasta_path: str,
    db_dir: str,
    max_kmer_len: int,
    charges: List[int],
) -> ProteinProductIonDb:
    # Set up
    t0 = time.time()
    logger.info(f"RUNNING create_protein_and_product_ion")
    # Check that files exist
    fasta_path = Path(fasta_path)
    db_path = Path(db_dir) / database_protein_product_ion_db_filename(
        fasta_path=fasta_path, max_k=max_kmer_len, charges_to_consider=charges
    )

    assert fasta_path.exists()
    if db_path.exists():
        print(f"Database path {db_path} already exists! Exiting and NOT overwriting!")
        sys.exit(0)
    else:
        try:
            create_protein_product_ion_db(
                fasta_path=fasta_path,
                db_path=db_path,
                max_kmer_len=max_kmer_len,
                charges_to_consider=charges,
            )
            print(f"Took {round(time.time() - t0, 2)} seconds")
        except Exception as e:
            # Delete database if it was created
            if db_path.exists():
                db_path.unlink()
            raise RuntimeError("Something failed!") from e


if __name__ == "__main__":
    create_protein_and_product_ion_database()
