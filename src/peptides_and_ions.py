import logging
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union

import click
from Bio import SeqIO
from fm_index import MultiFMIndex
from pydantic import BaseModel

from src.constants import (
    AMINO_ACID_MASSES,
    B_ION_TYPE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    HUMAN_PROTEOME,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
)
from src.utils import (
    ExistingPath,
    Kmer,
    PathType,
    Position,
    from_pickle,
    generate_aa_kmers,
    get_b_ion_prefixes,
    get_y_ion_suffixes,
    log_params,
    log_time,
    pickle_and_compress,
    setup_logger,
    to_pickle,
)

logger = logging.getLogger(__name__)


@dataclass
class ProteinRange(Position):
    protein: str

    @classmethod
    def from_pos(cls, protein: str, pos: Position):
        return cls(
            protein=protein,
            inclusive_start=pos.inclusive_start,
            exclusive_end=pos.exclusive_end,
        )

    def get_aa_seq(self, protein_name_to_seq_map: Dict[str, str]):
        return protein_name_to_seq_map[self.protein][
            self.inclusive_start : self.exclusive_end
        ]


class UnpositionedProductIon(BaseModel):
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

    @property
    def uid(self):
        return f"{self.ion_type}-{self.seq}-z{self.charge}"

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
                seq_generator = get_b_ion_prefixes
            elif ion_type == Y_ION_TYPE:
                seq_generator = get_y_ion_suffixes
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


class Peptide(BaseModel):
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
        self,
        charges: List[int],
        ion_types: Set[Literal[B_ION_TYPE, Y_ION_TYPE]] = {B_ION_TYPE, Y_ION_TYPE},
    ) -> List[UnpositionedProductIon]:
        return UnpositionedProductIon.generate_product_ions(
            seq=self.seq, charges=charges, ion_types=ion_types
        )

    def mz(self, charge: int) -> float:
        return compute_peptide_precursor_mz(seq=self.seq, charge=charge)


class Fasta(BaseModel):
    path: ExistingPath = HUMAN_PROTEOME

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

    def proteins_that_contain_seqs(self, seqs: List[str]) -> Dict[str, Set]:
        """Get the proteins that contain any of the query sequences."""
        seq_to_proteins = defaultdict(set)
        for protein in self.proteins:
            for seq in seqs:
                if seq in protein.seq:
                    seq_to_proteins[seq].add(protein.name)
        return dict(seq_to_proteins)

    @property
    def proteins(self) -> List[Peptide]:
        return Peptide.from_fasta(fasta_path=self.path)

    def get_proteins_by_name(self, names: Union[List[str], str, Path]) -> List[Peptide]:
        return get_proteins_by_name(protein_names=names, fasta_path=self.path)

    @staticmethod
    def write_fasta(peptides: List[Peptide], path: str) -> None:
        """
        Write a list of Peptide objects to a FASTA file.
        Each peptide will get two lines in the FASTA file:
        ><peptide.name> <peptide.desc>
        <peptide.seq>
        """
        with open(path, "w") as f:
            for peptide in peptides:
                header_parts = []
                if peptide.name:
                    header_parts.append(peptide.name)
                if peptide.desc:
                    header_parts.append(peptide.desc)

                header = " ".join(header_parts)
                f.write(f">{header}\n")
                f.write(f"{peptide.seq}\n")

    @cached_property
    def protein_name_to_seq_map(self):
        return {pep.name: pep.seq for pep in self.proteins}

    @cached_property
    def protein_name_to_peptide(self):
        return {pep.name: pep for pep in self.proteins}


@dataclass
class Fasta2MFMIndex:
    fasta: Fasta

    def __post_init__(self):
        if isinstance(self.fasta, (str, Path)):
            self.fasta = Fasta(path=self.fasta)

    def create_mfm_index(self) -> MultiFMIndex:
        fm_index = self.from_fasta(self.fasta.path)
        return fm_index

    @property
    def mfm_name(self) -> str:
        return self.fasta_name_to_mfm_index_name(fasta=self.fasta.path)

    @cached_property
    def idx_to_protein_map(self) -> Dict[int, Peptide]:
        return {idx: prot for idx, prot in enumerate(self.fasta.proteins)}

    @staticmethod
    def fasta_name_to_mfm_index_name(fasta: Union[str, Path]):
        return f"{Path(fasta).stem}.mfmindex"

    @staticmethod
    def mfm_index_name_to_fasta_name(mfm_index: Union[str, Path]):
        return f"{Path(mfm_index).stem}.fasta"

    @staticmethod
    def from_fasta(fasta_path: Union[str, Path]) -> MultiFMIndex:
        fasta = Fasta(path=fasta_path)
        seqs = [pep.seq for pep in fasta.proteins]
        return MultiFMIndex(data=seqs)

    @staticmethod
    def load(path: Union[str, Path]) -> MultiFMIndex:
        return from_pickle(path=path)

    @staticmethod
    def save(mfm_index: MultiFMIndex, path: Union[str, Path]):
        to_pickle(obj=mfm_index, path=path)

    @classmethod
    def create_and_save_index_from_fasta(
        cls, fasta: Union[str, Path], out_path: Union[str, Path]
    ):
        instance = cls(fasta=fasta)
        cls.save(mfm_index=instance.create_mfm_index(), path=out_path)


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
    pickle_and_compress(obj=uniq_kmers, path=output_path)


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
    name="create-mfm-index-for-fasta",
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
    "--overwrite",
    "-ow",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@click.option(
    "--out_dir",
    "-o",
    type=PathType(),
    required=False,
    help="If provided, the MFM index will be saved to <out_dir>/<FASTA file stem>.mfmindex. Otherwise it'll be saved to <FASTA file parent dir>/<FASTA file stem>.mfmindex.",
)
@log_time()
@log_params
def cli_create_mfm_index_for_fasta(
    fasta: Path,
    overwrite: bool,
    out_dir: Optional[Path],
):
    if out_dir is not None:
        out_path = out_dir / Fasta2MFMIndex.fasta_name_to_mfm_index_name(fasta=fasta)
    else:
        out_path = fasta.parent / Fasta2MFMIndex.fasta_name_to_mfm_index_name(
            fasta=fasta
        )
    if out_path.exists() and not overwrite:
        logger.info(
            f"File {out_path} already exists and overwrite={overwrite}. So skipping creation."
        )
    else:
        logger.info(f"Creating MFM index from FASTA and saving to {out_path}")
        mfm = Fasta2MFMIndex.from_fasta(fasta_path=fasta)
        Fasta2MFMIndex.save(mfm_index=mfm, path=out_path)


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


if __name__ == "__main__":
    logger = setup_logger()
    cli.add_command(cli_get_uniq_kmers)
    cli.add_command(cli_create_mfm_index_for_fasta)
    cli()
