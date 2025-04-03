from dataclasses import dataclass
from typing import Dict

from src.constants import AMINO_ACID_MASSES, MOUSE_PROTEOME
from src.peptides_and_ions import Peptide, get_uniq_kmer_to_protein_map
from src.sql_database import Sqlite3Database, SqlTableRow


@dataclass
class KmerWithMass(SqlTableRow):
    seq: str
    aa_mass: float

    @classmethod
    def from_seq(
        cls,
        seq: str,
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ):
        aa_mass = sum([amino_acid_mass_lookup[aa] for aa in seq])
        return cls(
            seq=seq,
            aa_mass=aa_mass,
        )


def test_smoke():
    # Constants
    num_prots = 2
    min_k, max_k = 1, 30
    table_name = "kmers"

    # Load proteins/amino acid sequences
    proteins = Peptide.from_fasta(MOUSE_PROTEOME)[:num_prots]

    # Get the unique kmers in the amino acid sequences and keep track of which
    # protein/amino acid sequence each kmer came from
    uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(
        proteins=proteins, min_k=min_k, max_k=max_k
    )

    # Form all hybrids

    # Create database from kmers
    db = Sqlite3Database()
    db.create_table_from_dataclass(table_name="kmers", obj=KmerWithMass)
    table_rows = []
    for kmer in uniq_kmer_to_protein_map.keys():
        table_rows.append(KmerWithMass.from_seq(seq=kmer))
    db.insert_dataclasses(table_name="kmers", data_classes=table_rows)
    db.add_index(table_name="kmers", index_name="mass", colms_to_index=["aa_mass"])
