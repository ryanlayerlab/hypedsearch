import logging
import pickle
import random
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from time import time
from typing import Dict, List, Optional, Set, Union

import click
from Bio import SeqIO
from pydantic import (BaseModel, field_validator)

from src.click_utils import PathType
from src.constants import (AMINO_ACID_MASSES, B_ION_TYPE, DEFAULT_MAX_K,
                           DEFAULT_MIN_K, PROTON_MASS, WATER_MASS,
                           Y_ION_TYPE, IonTypes)
from src.utils import (Kmer, generate_aa_kmers, get_time_in_diff_units,
                       setup_logger)

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

            header = " ".join(header_parts)
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


def compute_peptide_mz(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    return compute_y_ion_neutral_mass(
        aa_seq=aa_seq, charge=charge, amino_acid_mass_lookup=amino_acid_mass_lookup
    )


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
    """
    Get the proteins from a FASTA file.
    """
    proteins = []
    for p_id, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
        split_desc = protein.description.split(" ")
        name = split_desc[0]
        desc = " ".join(split_desc[1:])
        proteins.append(Peptide(seq=str(protein.seq), desc=desc, name=name, id=p_id))
    return proteins


@dataclass
class KmerToProteinIdMap:
    kmer_to_protein_id_map: Dict[str, List[int]]
    min_k: int = DEFAULT_MIN_K
    max_k: int = DEFAULT_MAX_K

    @classmethod
    def from_peptides(
        cls,
        peptides: List[Peptide],
        min_k: int = DEFAULT_MIN_K,
        max_k: int = DEFAULT_MAX_K,
    ):
        kmer_to_protein_map = get_uniq_kmer_to_protein_map(
            proteins=peptides, min_k=min_k, max_k=max_k
        )
        return cls(min_k=min_k, max_k=max_k, kmer_to_protein_id_map=kmer_to_protein_map)

    def save(self, path: Path):
        with open(path, "wb") as f:
            pickle.dump(
                self.kmer_to_protein_id_map, f, protocol=pickle.HIGHEST_PROTOCOL
            )

    @staticmethod
    def load(path: Path):
        with open(path, "rb") as f:
            return pickle.load(f)


class Fasta(BaseModel):
    fasta_path: Path

    @field_validator("fasta_path", mode="before")
    def str_to_path(cls, fasta_path: str) -> Path:
        return Path(fasta_path).absolute()

    @property
    def proteins(self) -> List[Peptide]:
        proteins = get_proteins_from_fasta(fasta_path=self.fasta_path)
        return proteins

    def proteins_by_name(self, protein_names: List[str]) -> List[Peptide]:
        proteins = get_proteins_by_name(
            protein_names=protein_names, fasta_path=self.fasta_path
        )
        return proteins


def get_proteins_by_name(
    protein_names: Union[List[str], Path],
    fasta_path: Optional[str] = None,
    proteins: Optional[List[Peptide]] = None,
) -> List[Peptide]:
    """
    Get proteins by name.
        - If fasta_path is provided, then the proteins returned are those in the fasta
        files with the given names (protein_names). protein_names can be specfied
        - protein_names can be provided as a list of names or as a Path to a file containing
        protein names
    """
    if fasta_path is not None:
        proteins = Peptide.from_fasta(fasta_path=fasta_path)

    if isinstance(protein_names, Path):
        protein_names = protein_names.read_text().splitlines()

    if protein_names is not None:
        proteins = list(filter(lambda protein: protein.name in protein_names, proteins))
        # assert len(proteins) == len(protein_names), (
        #     f"Expected there to be one row in the FASTA ({fasta_path}) "
        #     f"for each of the given proteins ({protein_names})"
        # )
    return proteins


def get_unique_kmers(peptides: Union[List[Peptide], List[str]], k: int) -> Set[str]:
    uniq_kmers = set()
    for peptide in peptides:
        if isinstance(peptide, str):
            peptide = Peptide(seq=peptide)
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
def get_uniq_kmer_to_protein_map(
    proteins: List[Peptide],
    min_k: int = DEFAULT_MIN_K,
    max_k: int = DEFAULT_MAX_K,
) -> Dict[str, List[int]]:
    """ """
    fcn_start_time = time()
    uniq_kmer_to_protein_map = defaultdict(list)
    num_proteins = len(proteins)
    for p_idx, protein in enumerate(proteins):
        if p_idx % 10 == 0:
            logger.info(f"Processing protein {p_idx+1} of {num_proteins}")
        prot_start_time = time()
        uniq_kmers = set(kmer.seq for kmer in protein.kmers(min_k=min_k, max_k=max_k))
        for kmer in uniq_kmers:
            uniq_kmer_to_protein_map[kmer].append(protein.id)
        # logger.info(f"\t took {round(time()-prot_start_time, 2)} seconds")
    total_time = time() - fcn_start_time
    logger.info(
        f"Creating the unique kmers-to-protein map took {get_time_in_diff_units(time_sec=total_time)}"
    )
    logger.info(f"Number of unique kmers {len(uniq_kmer_to_protein_map)}")
    return uniq_kmer_to_protein_map


def create_kmer_to_protein_id_map(
    fasta_path: Path,
    min_k: int,
    max_k: int,
    output_path: Optional[Path] = None,
    protein_names: Optional[Path] = None,
):
    logger.info("Creating kmer-to-protein-id map...")
    start_time = time()

    # Load proteins
    proteins = Peptide.from_fasta(fasta_path=fasta_path)
    if protein_names is not None:
        proteins = get_proteins_by_name(proteins=proteins, protein_names=protein_names)

    # Create kmer-to-protein-id map
    k_to_p_map = KmerToProteinIdMap.from_peptides(
        peptides=proteins, min_k=min_k, max_k=max_k
    )

    # Return map if you don't want to save it
    if output_path is None:
        return k_to_p_map

    # Save map
    if output_path.is_dir():
        if protein_names is not None:
            output_path = (
                output_path
                / f"{fasta_path.stem}.{protein_names.stem}.kmer_to_protein_map.pkl"
            )
        else:
            output_path = output_path / f"{fasta_path.stem}.kmer_to_protein_map.pkl"

    k_to_p_map.save(path=output_path)

    logger.info(
        f"Done creating kmer-to-protein-id map. It took {get_time_in_diff_units(time() - start_time)}"
    )
    return k_to_p_map


@click.command(
    name="create-kmer-to-protein-id-map",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--fasta_path",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--output_path",
    "-o",
    type=PathType(),
    required=True,
    help=(
        "Where to save the kmer-to-protein-id map. "
        "If a directory is provided, the map will be saved as "
        "<FASTA file stem>.<protein_names stem if provided>.kmer_to_protein_map.pkl. "
        "If a path to a .pkl file is provided, the map will be saved to that path."
    ),
)
@click.option(
    "--min_k",
    "-mk",
    type=int,
    default=DEFAULT_MIN_K,
    show_default=True,
    help="Minimum kmer length to consider.",
)
@click.option(
    "--max_k",
    "-Mk",
    type=int,
    default=DEFAULT_MAX_K,
    show_default=True,
    help="Maximum kmer length to consider.",
)
@click.option(
    "--protein_names",
    "-pn",
    type=PathType(),
    required=False,
    help=(
        "Path to file containing the names of proteins to form hybrids from. "
        "If not provided, all proteins in the FASTA will be used"
    ),
)
def cli_create_kmer_to_protein_id_map(
    fasta_path: Path,
    min_k: int,
    max_k: int,
    output_path: Path,
    protein_names: Optional[Path] = None,
):
    create_kmer_to_protein_id_map(
        fasta_path=fasta_path,
        min_k=min_k,
        max_k=max_k,
        output_path=output_path,
        protein_names=protein_names,
    )


if __name__ == "__main__":
    setup_logger()
    cli_create_kmer_to_protein_id_map()
