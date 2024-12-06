import logging
import sqlite3
import time
from dataclasses import asdict
from typing import Dict, Iterator, List

import click

from src.erik import generate_kmers_with_masses
from src.erik_constants import CHARGE, LOGS_DIR, MASS, MAX_PEPTIDE_LEN, SEQ
from src.erik_utils import Protein, get_proteins_from_fasta, log_params, setup_logger
from src.lookups.constants import AMINO_ACID_MASSES

setup_logger(log_file=LOGS_DIR / f"{__name__}.log")
logger = logging.getLogger(__name__)


def create_kmer_mass_db(
    proteins: Iterator[Protein],
    db_file: str,
    table_name: str = "kmers",
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    max_peptide_len: float = MAX_PEPTIDE_LEN,
    charges: List[int] = [1],
) -> None:

    # Create DB
    connection = sqlite3.connect(db_file, timeout=10)
    cursor = connection.cursor()

    _ = cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
    _ = cursor.execute(
        f"""
        CREATE TABLE {table_name} (
            {SEQ} TEXT PRIMARY KEY,
            {MASS} REAL,
            {CHARGE} INTEGER 
        )
        """
    )
    t0 = time.time()
    for p_num, protein in enumerate(proteins):
        # Log progress
        if (p_num + 1) % 10 == 0:
            logger.info(f"Protein {p_num+1}; {round(time.time()-t0, 2)} seconds")

        # Generate kmers from protein and add them to DB if they don't already exist
        kmers = generate_kmers_with_masses(
            peptide=protein.seq,
            max_kmer_len=max_peptide_len,
            amino_acid_mass_lookup=amino_acid_mass_lookup,
        )
        # Incorporate charge
        kmers = [
            (kmer.seq, kmer.mass * charge, charge)
            for kmer in kmers
            for charge in charges
        ]
        _ = cursor.executemany(
            f"INSERT OR IGNORE INTO {table_name} VALUES(?, ?, ?)",
            kmers,
        )
        _ = connection.commit()

    connection.close()


@click.command(
    help="""
    This script creates a database with a table with colums:\n
    <kmer seq>, <kmer mass>, <charge>\n
    """,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--fasta_path", "-f", type=str, required=True, help="Path to FASTA file")
@click.option(
    "--db_file",
    "-d",
    type=str,
    required=True,
    help="Path to to-be-created kmer database",
)
@click.option(
    "--max_kmer_len",
    "-k",
    help=f"Maximum kmer length to add to database. Defaults to {MAX_PEPTIDE_LEN}",
    type=int,
    required=False,
    default=MAX_PEPTIDE_LEN,
)
@click.option(
    "--charges",
    "-c",
    type=int,
    multiple=True,
    required=False,
    help="A list of integer charges, e.g., --charges 1 2 3. Defaults to [1]",
)
@log_params
def main(
    fasta_path: str,
    db_file: str,
    table_name: str = "kmers",
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    max_kmer_len: float = MAX_PEPTIDE_LEN,
    charges: List[int] = [1],
) -> None:
    try:
        logger.info(f"Getting proteins...")
        proteins = get_proteins_from_fasta(fasta_path=fasta_path)
        logger.info("Creating kmer database...")
        create_kmer_mass_db(
            proteins=list(proteins),
            db_file=db_file,
            table_name=table_name,
            amino_acid_mass_lookup=amino_acid_mass_lookup,
            max_peptide_len=max_kmer_len,
            charges=charges,
        )
        logger.info("Finshed creating db. Exiting")
    except:
        logger.info("An error happened somewhere! Exiting")


if __name__ == "__main__":
    main()
