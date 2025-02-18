from typing import Generator, List

from Bio import SeqIO

from src.lookups.data_classes import Protein


def get_proteins_from_fasta(fasta_path: str) -> Generator[Protein, None, None]:
    for protein_num, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
        yield Protein(
            protein_id=protein_num,
            sequence=str(protein.seq),
            desc=protein.description,
        )


def get_specific_protein_from_fasta(fasta_path: str, protein_name: str) -> Protein:
    proteins = list(get_proteins_from_fasta(fasta_path=fasta_path))
    matching_proteins = list(
        filter(lambda protein: protein.desc.split(" ")[0] == protein_name, proteins)
    )
    assert len(matching_proteins) == 1, (
        f"Expected there to be one row in the FASTA ({fasta_path}) "
        f"with the given protein name ({protein_name})"
    )
    protein = matching_proteins[0]
    return protein
