from typing import Generator

from Bio import SeqIO

from src.lookups.data_classes import Protein


def get_proteins_from_fasta(fasta_path: str) -> Generator[Protein, None, None]:
    for protein_num, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
        yield Protein(
            protein_id=protein_num,
            sequence=str(protein.seq),
            desc=protein.description,
        )
