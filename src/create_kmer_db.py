import logging
import sqlite3
import time
from dataclasses import asdict
from typing import Dict, Iterator

import click

from src.erik import generate_kmers_with_masses
from src.erik_constants import LOGS_DIR, MASS, MAX_PEPTIDE_LEN, SEQ
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
) -> None:

    # Create DB
    connection = sqlite3.connect(db_file, timeout=10)
    cursor = connection.cursor()

    _ = cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
    _ = cursor.execute(
        f"""
        CREATE TABLE {table_name} (
            {SEQ} TEXT PRIMARY KEY,
            {MASS} REAL
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
        _ = cursor.executemany(
            f"INSERT OR IGNORE INTO {table_name} VALUES(:seq, :mass)",
            [asdict(kmer) for kmer in list(kmers)],
        )
        _ = connection.commit()

    connection.close()


@click.command(help="hello!", context_settings=dict(help_option_names=["-h", "--help"]))
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
    help="Maximum kmer length to add to database",
    type=int,
    required=False,
    default=MAX_PEPTIDE_LEN,
)
@log_params
def main(
    fasta_path: str,
    db_file: str,
    table_name: str = "kmers",
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    max_kmer_len: float = MAX_PEPTIDE_LEN,
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
        )
        logger.info("Finshed creating db. Exiting")
    except:
        logger.info("An error happened somewhere! Exiting")


if __name__ == "__main__":
    main()
