import logging
import random
from abc import ABC, abstractmethod
from collections import defaultdict
from time import time
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Callable, Dict, List, Literal, Optional, Set

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
    IonTypes,
)
from src.utils import Kmer, Position, generate_aa_kmers, get_time_in_diff_units, log_params, to_path

logger = logging.getLogger(__name__)


@dataclass
class ProductIon:
    seq: str
    charge: int
    ion_type: IonTypes
    neutral_mass: float = field(init=False)

    def __post_init__(self):
        if self.ion_type == IonTypes.B_ION_TYPE:
            self.neutral_mass = compute_b_ion_neutral_mass(
                aa_seq=self.seq, charge=self.charge
            )
        elif self.ion_type == IonTypes.Y_ION_TYPE:
            self.neutral_mass = compute_y_ion_neutral_mass(
                aa_seq=self.seq, charge=self.charge
            )
        else:
            msg = (
                "Can't compute ion neutral mass because the ion type isn't supported. "
                f"Received ion type {self.ion_type}"
            )
            raise RuntimeError(msg)

        self.ion_type_as_str = self.ion_type.value


@dataclass
class Peptide:
    seq: str
    name: Optional[str] = None
    desc: Optional[str] = None
    id: Optional[int] = None

    @classmethod
    def from_fasta(cls, fasta_path: str) -> List["Peptide"]:
        return get_proteins_from_fasta(fasta_path=fasta_path)

    def kmers(self, max_k: int, min_k: int = 1) -> List[Kmer]:
        return generate_aa_kmers(aa_seq=self.seq, min_k=min_k, max_k=max_k)

    def product_ions(
        self, ion_types: List[IonTypes], charges: List[int]
    ) -> List[ProductIon]:
        return generate_product_ions(seq=self.seq, charges=charges, ion_types=ion_types)


def write_fasta(peptides: List[Peptide], output_path: str) -> None:
    with open(output_path, "w") as f:
        for peptide in peptides:
            header_parts = []
            if peptide.name:
                header_parts.append(peptide.name)
            # if peptide.id is not None:
            #     header_parts.append(f"id:{peptide.id}")
            if peptide.desc:
                header_parts.append(peptide.desc)

            header = " ".join(header_parts) or "unnamed_peptide"
            f.write(f">{header}\n")
            f.write(f"{peptide.seq}\n")


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


def get_product_ion_creator(ion_type: IonTypes) -> ProductIonCreator:
    if ion_type == IonTypes.B_ION_TYPE:
        return BIonCreator
    elif ion_type == IonTypes.Y_ION_TYPE:
        return YIonCreator
    else:
        msg = (
            "Can't compute ion neutral mass because the ion type isn't supported.\n"
            f"Supported types are [{B_ION_TYPE}, {Y_ION_TYPE}].\n"
            f"Received ion type {ion_type}"
        )
        raise RuntimeError(msg)


def generate_product_ions(seq: str, charges: List[int], ion_types: List[IonTypes]):
    product_ions = []
    for ion_type in ion_types:
        seq_generator = get_product_ion_creator(
            ion_type=ion_type
        ).generate_product_ion_seqs
        for charge in charges:
            product_ions.extend(
                [
                    ProductIon(seq=ion_seq, charge=charge, ion_type=ion_type)
                    for ion_seq in seq_generator(seq=seq)
                ]
            )
    return product_ions


def get_proteins_from_fasta(fasta_path: str) -> List[Peptide]:
    proteins = []
    for p_id, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
        desc = protein.description
        name = desc.split(" ")[0]
        proteins.append(Peptide(seq=str(protein.seq), desc=desc, name=name, id=p_id))
    return proteins


def get_proteins_by_name(
    protein_names: str,
    fasta_path: Optional[str] = None,
    proteins: Optional[List[Peptide]] = None,
) -> List[Peptide]:
    if fasta_path is not None:
        proteins = list(get_proteins_from_fasta(fasta_path=fasta_path))
    matching_proteins = list(
        filter(lambda protein: protein.desc.split(" ")[0] in protein_names, proteins)
    )
    assert len(matching_proteins) == len(protein_names), (
        f"Expected there to be one row in the FASTA ({fasta_path}) "
        f"for each of the given proteins ({protein_names})"
    )
    return matching_proteins


def get_unique_kmers(peptides: List[Peptide], k: int) -> Set[str]:
    uniq_kmers = set()
    for peptide in peptides:
        peptide_kmers = {kmer.seq for kmer in peptide.kmers(min_k=k, max_k=k)}
        uniq_kmers.update(peptide_kmers)
    return uniq_kmers


def random_sample_of_unique_kmers(
    k: int,
    sample_size: int,
    peptides: Optional[List[Peptide]] = None,
    fasta_path: Optional[str] = None,
) -> List[str]:
    if fasta_path is not None:
        peptides = Peptide.from_fasta(fasta_path=fasta_path)

    uniq_kmers = get_unique_kmers(peptides=peptides, k=k)
    if len(uniq_kmers) < sample_size:
        raise RuntimeError(
            f"The number of unique kmers, {len(uniq_kmers)}, must be >= than the sample size, {sample_size}"
        )

    return random.sample(sorted(uniq_kmers), k=sample_size)

# Don't use @log_params because 'proteins' can be a list of many, many peptides
def get_unique_peptides(
    min_k: int,
    max_k: int,
    proteins: List[Peptide],
):
    fcn_start_time = time()
    uniq_peptides = defaultdict(list)
    num_proteins = len(proteins)
    for p_idx, protein in enumerate(proteins):
        logger.info(f"Processing protein {p_idx+1} of {num_proteins}")
        prot_start_time = time()
        uniq_kmers = set(kmer.seq for kmer in protein.kmers(min_k=min_k, max_k=max_k))
        for kmer in uniq_kmers:
            uniq_peptides[kmer].append(protein.id)
        logger.info(f"\t took {round(time()-prot_start_time, 2)} seconds")
    total_time = time() - fcn_start_time
    logger.info(f"In total, function took {get_time_in_diff_units(time_sec=total_time)}")
    return uniq_peptides
