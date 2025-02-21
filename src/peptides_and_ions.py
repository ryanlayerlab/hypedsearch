from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Literal, Optional

from Bio import SeqIO
from pydantic import BaseModel, BeforeValidator, computed_field
from typing_extensions import Annotated

from src.constants import (
    ALL_IONS,
    AMINO_ACID_MASSES,
    B_ION_TYPE,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
)
from src.utils import Kmer, Position, generate_aa_kmers, to_path


class Fasta(BaseModel):
    path: Annotated[Path, BeforeValidator(to_path)]

    def proteins(self):
        return get_proteins_from_fasta(fasta_path=self.path)


class ProductIon(BaseModel):
    seq: str
    charge: int
    ion_type: Literal[B_ION_TYPE, Y_ION_TYPE]

    @computed_field
    @property
    def neutral_mass(self) -> float:
        if self.ion_type == B_ION_TYPE:
            return compute_b_ion_neutral_mass(aa_seq=self.seq, charge=self.charge)
        elif self.ion_type == Y_ION_TYPE:
            return compute_y_ion_neutral_mass(aa_seq=self.seq, charge=self.charge)
        else:
            msg = (
                "Can't compute ion neutral mass because the ion type isn't supported. "
                f"Received ion type {self.ion_type}"
            )
            raise RuntimeError(msg)


class Peptide(BaseModel):
    seq: str
    name: Optional[str] = None
    desc: Optional[str] = None
    id: Optional[int] = None

    @classmethod
    def from_fasta(cls, fasta_path: str) -> List[Peptide]:
        fasta = Fasta(path=fasta_path)
        return fasta.proteins()

    def kmers(self, max_k: int, min_k: Optional[int] = 1) -> List[Kmer]:
        return generate_aa_kmers(aa_seq=self.seq, min_k=min_k, max_k=max_k)

    def product_ions(
        self, ion_type: Literal[B_ION_TYPE, Y_ION_TYPE, ALL_IONS], charges: List[int]
    ) -> List[ProductIon]:
        return generate_product_ions(seq=self.seq, charges=charges, ion_type=ion_type)


def compute_b_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + (charge * PROTON_MASS)) / charge
    return neutral_mass


def compute_y_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + WATER_MASS + (charge * PROTON_MASS)) / charge
    return neutral_mass


class ProductIonCreator(ABC):
    @staticmethod
    @abstractmethod
    def generate_product_ion_seqs(seq: str) -> List[str]:
        pass


class BIonCreator(ProductIonCreator):
    @staticmethod
    def generate_product_ion_seqs(seq: str) -> List[str]:
        return [seq[:i] for i in range(1, len(seq) + 1)]  # Prefixes (b-ions)


class YIonCreator(ProductIonCreator):
    @staticmethod
    def generate_product_ion_seqs(seq: str) -> List[str]:
        return [seq[i:] for i in range(0, len(seq))]  # Suffixes (y-ions)


def get_product_ion_seq_generator(ion_type: str):
    if ion_type == B_ION_TYPE:
        return BIonCreator.generate_product_ion_seqs
    elif ion_type == Y_ION_TYPE:
        return YIonCreator.generate_product_ion_seqs
    else:
        msg = (
            "Can't compute ion neutral mass because the ion type isn't supported.\n"
            f"Supported types are [{B_ION_TYPE}, {Y_ION_TYPE}].\n"
            f"Received ion type {ion_type}"
        )
        raise RuntimeError(msg)


def generate_product_ions(seq: str, charges: List[int], ion_type: str):
    seq_generator = get_product_ion_seq_generator(ion_type=ion_type)

    product_ions = []
    for charge in charges:
        product_ions.extend(
            [
                ProductIon(seq=ion_seq, charge=charge, ion_type=ion_type)
                for ion_seq in seq_generator(seq=seq)
            ]
        )
    return product_ions


def get_proteins_from_fasta(fasta_path: str) -> List[Peptide]:
    proteins = [
        Peptide(seq=str(protein.seq), desc=protein.description, id=p_id)
        for p_id, protein in enumerate(SeqIO.parse(fasta_path, "fasta"))
    ]
    return proteins


def get_specific_protein_from_fasta(fasta_path: str, protein_name: str) -> Peptide:
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
