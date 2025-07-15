import logging
from pathlib import Path
from time import time
from typing import List, Optional, Union

import click

from src.constants import (
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    MEMORY,
    MOUSE_PROTEOME,
)
from src.peptides_and_ions import (
    Peptide,
    get_proteins_by_name,
    get_uniq_kmer_to_protein_map,
)
from src.protein_product_ion_database import DbKmer, DbProtein, ProteinProductIonDb
from src.utils import PathType, get_time_in_diff_units, setup_logger

logger = logging.getLogger(__name__)


def get_proteins_for_db(
    proteins: Optional[Union[List[Peptide], Path, List[str]]], fasta_path: Path
) -> List[Peptide]:
    # Case 1: proteins is Path to a file containing protein names
    if isinstance(proteins, Path):
        with open(proteins, "r") as f:
            protein_names = [line.strip() for line in f.readlines()]
        db_proteins = get_proteins_by_name(
            protein_names=protein_names, fasta_path=fasta_path
        )
    # Case 2: proteins is a list
    elif isinstance(proteins, list):
        # Case 2a: proteins are Peptides
        if isinstance(proteins[0], Peptide):
            db_proteins = proteins
        # Case 2b: proteins are names of proteins
        elif isinstance(proteins[0], str):
            db_proteins = get_proteins_by_name(
                protein_names=proteins, fasta_path=fasta_path
            )
    # Case 3: proteins is None so get all proteins from the FASTA file
    elif proteins is None:
        # Assuming fasta_path is defined globally or passed in some way
        db_proteins = Peptide.from_fasta(fasta_path=fasta_path)

    # Case 4: proteins is not a valid type
    else:
        raise ValueError(
            "Proteins must be either a Path to a file containing protein names, "
            "a list of Peptide objects, or None."
        )
    return db_proteins


def create_protein_table(db: ProteinProductIonDb, db_proteins: List[Peptide]) -> None:
    db_proteins = [DbProtein.from_peptide(prot) for prot in db_proteins]
    db.create_table_from_dataclass(table_name=db.protein_table_name, obj=DbProtein)
    db.insert_dataclasses(table_name=db.protein_table_name, data_classes=db_proteins)


def create_product_ion_table(
    db: ProteinProductIonDb,
    db_proteins: List[Peptide],
    min_k: int = DEFAULT_MIN_KMER_LEN,
    max_k: int = DEFAULT_MAX_KMER_LEN,
) -> None:
    uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(
        proteins=db_proteins, min_k=min_k, max_k=max_k
    )
    logger.info("Creating 'DbKmer' objects for the database")
    ions = [
        DbKmer.seq_to_ion(seq=seq, protein_ids=protein_ids)
        for seq, protein_ids in uniq_kmer_to_protein_map.items()
    ]
    logger.info(
        "Creating the table for the 'DbKmer' object and inserting the object instances"
    )
    db.create_table_from_dataclass(table_name=db.product_ion_table_name, obj=DbKmer)
    db.insert_dataclasses(table_name=db.product_ion_table_name, data_classes=ions)
    logger.info("Creating 'mass' index for the table")
    db.add_index(
        table_name=db.product_ion_table_name,
        index_name="mass",
        colms_to_index=["aa_mass"],
    )


def create_db(
    fasta_path: Optional[Path] = None,
    proteins: Optional[Union[List[Peptide], Path]] = None,
    db_path: Union[Path, str] = MEMORY,
    min_k: int = DEFAULT_MIN_KMER_LEN,
    max_k: int = DEFAULT_MAX_KMER_LEN,
) -> ProteinProductIonDb:
    t0 = time()
    logger.info(f"Creating DB at {db_path}")
    # Initialize the database
    db = ProteinProductIonDb(
        db_path=db_path,
    )

    # Create protein table
    db_proteins = get_proteins_for_db(proteins=proteins, fasta_path=fasta_path)
    create_protein_table(db=db, db_proteins=db_proteins)

    # Create product/fragment ion table
    create_product_ion_table(
        db=db,
        db_proteins=db_proteins,
        min_k=min_k,
        max_k=max_k,
    )

    logger.info(
        f"Creating and indexing DB with {len(db_proteins)} took {get_time_in_diff_units(time() - t0)}"
    )

    return db


@click.command(
    name="create-db",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="Create a database of product ions",
)
@click.option(
    "--db_path",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the to-be-created database.",
)
@click.option(
    "--fasta_path",
    "-f",
    type=PathType(),
    default=MOUSE_PROTEOME,
    show_default=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--proteins",
    "-p",
    type=PathType(),
    required=False,
    help=(
        "Path to the file containing protein names, one per line. "
        "The protein names should be in the FASTA file. "
        "If not provided, all proteins from the FASTA file will be used."
    ),
)
@click.option(
    "--min_k",
    "-mk",
    type=int,
    default=DEFAULT_MIN_KMER_LEN,
    show_default=True,
    help="Minimum kmer length to consider.",
)
@click.option(
    "--max_k",
    "-Mk",
    type=int,
    default=DEFAULT_MAX_KMER_LEN,
    show_default=True,
    help="Maximum kmer length to consider.",
)
def create_db_cli(
    fasta_path: Path,
    db_path: Path,
    proteins: Optional[Union[List[Peptide], Path]] = None,
    min_k: int = DEFAULT_MIN_KMER_LEN,
    max_k: int = DEFAULT_MAX_KMER_LEN,
):
    create_db(
        fasta_path=fasta_path,
        proteins=proteins,
        db_path=db_path,
        min_k=min_k,
        max_k=max_k,
    )


if __name__ == "__main__":
    setup_logger()
    create_db_cli()
