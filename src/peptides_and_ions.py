import logging
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from time import time
from typing import Dict, List, Literal, Optional, Set, Union

import click
from Bio import SeqIO
from pydantic import BaseModel

from src.constants import (
    AMINO_ACID_MASSES,
    B_ION_TYPE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
    IonTypes,
)
from src.utils import (
    ExistingPath,
    Kmer,
    PathType,
    generate_aa_kmers,
    get_time_in_diff_units,
    log_params,
    log_time,
    pickle_and_compress,
    prefixes,
    setup_logger,
    suffixes,
)

logger = logging.getLogger(__name__)


@dataclass
class UnpositionedProductIon:
    seq: str
    charge: int
    ion_type: Literal[B_ION_TYPE, Y_ION_TYPE]
    proteins: Optional[List[Union[str, int]]] = None

    @property
    def mz(self) -> float:
        return self.compute_ion_mz(
            seq=self.seq,
            charge=self.charge,
            ion_type=self.ion_type,
        )

    @staticmethod
    def compute_b_ion_mz(
        seq: str,
        charge: int,
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ) -> float:
        aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in seq])
        neutral_mass = (aa_mass_sum + (charge * PROTON_MASS)) / charge
        return neutral_mass

    @staticmethod
    def compute_y_ion_mz(
        seq: str,
        charge: int,
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ) -> float:
        aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in seq])
        neutral_mass = (aa_mass_sum + WATER_MASS + (charge * PROTON_MASS)) / charge
        return neutral_mass

    @staticmethod
    def compute_ion_mz(
        seq: str,
        charge: int,
        ion_type: Literal[B_ION_TYPE, Y_ION_TYPE],
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ):
        if ion_type == B_ION_TYPE:
            return UnpositionedProductIon.compute_b_ion_mz(
                seq=seq,
                charge=charge,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
            )
        elif ion_type == Y_ION_TYPE:
            return UnpositionedProductIon.compute_y_ion_mz(
                seq=seq,
                charge=charge,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
            )
        else:
            raise ValueError(f"Unsupported ion type: {ion_type}")

    @classmethod
    def generate_product_ions(
        cls,
        seq: str,
        charges: List[int],
        ion_types: Set[Literal[B_ION_TYPE, Y_ION_TYPE]] = {B_ION_TYPE, Y_ION_TYPE},
    ) -> List["UnpositionedProductIon"]:
        product_ions = []
        for ion_type in ion_types:
            if ion_type == B_ION_TYPE:
                seq_generator = prefixes
            elif ion_type == Y_ION_TYPE:
                seq_generator = suffixes
            else:
                raise ValueError(f"Unsupported ion type: {ion_type}")
            for charge in charges:
                product_ions.extend(
                    [
                        cls(seq=ion_seq, charge=charge, ion_type=ion_type)
                        for ion_seq in seq_generator(seq)
                    ]
                )
        return product_ions


@dataclass
class Peptide:
    seq: str
    name: Optional[str] = None
    desc: Optional[str] = None
    id: Optional[int] = None

    @classmethod
    def from_fasta(cls, fasta_path: str) -> List["Peptide"]:
        proteins = []
        for p_id, protein in enumerate(SeqIO.parse(fasta_path, "fasta")):
            split_desc = protein.description.split(" ")
            name = split_desc[0]
            desc = " ".join(split_desc[1:])
            proteins.append(cls(seq=str(protein.seq), desc=desc, name=name, id=p_id))
        return proteins

    def kmers(self, max_k: int, min_k: int = 1) -> List[Kmer]:
        return generate_aa_kmers(aa_seq=self.seq, min_k=min_k, max_k=max_k)

    def product_ions(
        self, charges: List[int], ion_types: List[IonTypes] = {B_ION_TYPE, Y_ION_TYPE}
    ) -> List[UnpositionedProductIon]:
        return UnpositionedProductIon.generate_product_ions(
            seq=self.seq, charges=charges, ion_types=ion_types
        )


class Fasta(BaseModel):
    path: ExistingPath

    @cached_property
    def seqs(self):
        """Get the sequences from the FASTA file."""
        return [str(record.seq) for record in SeqIO.parse(self.path, "fasta")]

    def contains_seq(self, query_seq: str) -> bool:
        """Check if the query sequence exists in the FASTA file."""
        for seq in self.seqs:
            if query_seq in seq:
                return True
        return False

    @property
    def proteins(self) -> List[Peptide]:
        return Peptide.from_fasta(fasta_path=self.path)

    def get_proteins_by_name(
        self, protein_names: Union[List[str], str, Path]
    ) -> List[Peptide]:
        return get_proteins_by_name(protein_names=protein_names, fasta_path=self.path)

    @staticmethod
    def write_fasta(peptides: List[Peptide], out_path: str) -> None:
        """
        Write a list of Peptide objects to a FASTA file.
        Each peptide will get two lines in the FASTA file:
        ><peptide.name> <peptide.desc>
        <peptide.seq>
        """
        with open(out_path, "w") as f:
            for peptide in peptides:
                header_parts = []
                if peptide.name:
                    header_parts.append(peptide.name)
                if peptide.desc:
                    header_parts.append(peptide.desc)

                header = " ".join(header_parts)
                f.write(f">{header}\n")
                f.write(f"{peptide.seq}\n")

    @property
    def protein_name_to_seq_map(self):
        return {pep.name: pep.seq for pep in self.proteins}


def compute_peptide_precursor_mz(seq: str, charge: int):
    """
    The m/z of a peptide as a precursor
        peptide_mz = (sum_AA + WATER + z*PROTON) / z
    which is the same as the same as if the sequence is considered a y-ion
    """
    return UnpositionedProductIon.compute_y_ion_mz(seq=seq, charge=charge)


def get_proteins_by_name(
    protein_names: Union[List[str], Path],
    fasta_path: Optional[str] = None,
    proteins: Optional[List[Peptide]] = None,
) -> List[Peptide]:
    if fasta_path is not None:
        proteins = Peptide.from_fasta(fasta_path=fasta_path)
    if isinstance(protein_names, Path):
        protein_names = protein_names.read_text().splitlines()
    if protein_names is not None:
        proteins = list(filter(lambda protein: protein.name in protein_names, proteins))
    return proteins


def get_unique_kmers(
    peptides: Union[List[Peptide], List[str], Path], min_k: int, max_k: int
) -> Set[str]:
    """
    Given a list of amino acid sequences (either as strings, Peptide objects, or a path to a FASTA file),
    return all unique k-mers from the sequences for k=min_k, min_k + 1, ..., max_k
    """
    if isinstance(peptides, Path):
        logger.info("Reading in FASTA file...")
        peptides = Peptide.from_fasta(fasta_path=peptides)
        logger.info("Done reading in FASTA file")
    uniq_kmers = set()
    num_proteins = len(peptides)
    logger.info(f"Number of proteins: {num_proteins}")
    for p_idx, peptide in enumerate(peptides):
        if p_idx % 100 == 0:
            logger.info(f"Processing protein {p_idx + 1} of {num_proteins}")
        if isinstance(peptide, str):
            peptide = Peptide(seq=peptide)
        uniq_kmers.update(
            {kmer.seq for kmer in peptide.kmers(min_k=min_k, max_k=max_k)}
        )
    return uniq_kmers


@click.command(
    name="get-uniq-kmers",
    help="Get all unique kmers from a FASTA file. Resulting set object will be pickled and compresses",
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
@log_time(level=logging.DEBUG)
@log_params
def cli_get_uniq_kmers(
    fasta_path: Path,
    min_k: int,
    max_k: int,
    output_path: Path,
):
    """
    Get all unique kmers from a FASTA file. Resulting set object will be pickled and compresses
    """
    logger.info("Getting unique kmers...")
    uniq_kmers = get_unique_kmers(peptides=fasta_path, min_k=min_k, max_k=max_k)
    logger.info(f"Number of unique kmers: {len(uniq_kmers)}")

    logger.info("Pickling and compressing unique kmers...")
    if output_path.is_dir():
        output_path = output_path / f"{fasta_path.stem}.uniq_kmers.pklz"
    pickle_and_compress(obj=uniq_kmers, file_path=output_path)


def get_kmer_counts_by_protein(
    fasta: Path,
    k: int,
) -> Dict[str, Dict[str, int]]:
    """
    Given a FASTA file and a k, return a dictionary that looks like the following:
    {<kmer sequence>: {<protein name>: <number of time kmer appears in protein>}}
    """
    peptides = Peptide.from_fasta(fasta_path=fasta)
    kmer_to_prot_to_count_map = defaultdict(lambda: defaultdict(int))
    num_prots = len(peptides)
    for idx, peptide in enumerate(peptides):
        logger.debug(f"Processing protein {idx + 1} of {num_prots}")
        kmers = peptide.kmers(min_k=k, max_k=k)
        num_kmers = len(kmers)
        logger.debug(f"Number of kmers: {num_kmers}")
        for kmer in kmers:
            kmer_to_prot_to_count_map[kmer.seq][peptide.name] += 1

    return dict(kmer_to_prot_to_count_map)


@click.command(
    name="get-kmer-info",
    help=(
        "Given a FASTA and a k, get the number of times each unique k-mer appears in each protein."
    ),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--k",
    "-k",
    type=int,
    required=True,
    help="The k-mer length, k",
)
@click.option(
    "--out_path",
    "-o",
    type=PathType(),
    required=True,
    help=("The results will be saved in this location"),
)
@click.option(
    "--overwrite",
    "-ow",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@log_time(level=logging.DEBUG)
@log_params
def cli_get_kmer_counts_by_protein(
    fasta: Path,
    k: int,
    out_path: Path,
    overwrite: bool,
):
    if out_path.exists() and not overwrite:
        logger.info(f"File {out_path} already exists. Skipping...")
        return
    loop_start_time = time()
    logger.info(f"Processing k={k}")
    k_map = get_kmer_counts_by_protein(
        fasta=fasta,
        k=k,
    )
    pickle_and_compress(obj=k_map, file_path=out_path)
    loop_duration = time() - loop_start_time
    logger.info(f"Processed k={k} in {get_time_in_diff_units(loop_duration)}")


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


if __name__ == "__main__":
    logger = setup_logger()
    cli.add_command(cli_get_uniq_kmers)
    cli.add_command(cli_get_kmer_counts_by_protein)
    cli()
