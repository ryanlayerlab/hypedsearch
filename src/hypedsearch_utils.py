import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

from src.constants import HS_PREFIX
from src.peptides_and_ions import (
    Fasta,
    Peptide,
    compute_peptide_mz,
    get_proteins_by_name,
    get_uniq_kmer_to_protein_map,
    write_fasta,
)
from src.sql_database import SqlTableRow
from src.utils import load_json, log_time

logger = logging.getLogger(__name__)


@dataclass
class HybridPeptide:
    b_seq: str
    y_seq: str
    b_prot_ids: Set[int] = field(default_factory=set)
    y_prot_ids: Set[int] = field(default_factory=set)
    b_prot_names: Optional[Set[str]] = None
    y_prot_names: Optional[Set[str]] = None
    fasta_description: Optional[str] = field(init=False, default=None)
    scan: Optional[int] = None
    sample: Optional[str] = None

    def __post_init__(self):
        if (self.b_prot_names is not None) and (self.y_prot_names is not None):
            self.fasta_description = f"b-prots:{','.join(self.b_prot_names)} y-prots:{','.join(self.y_prot_names)}"

    @property
    def seq(self):
        return self.b_seq + self.y_seq

    @property
    def seq_with_hyphen(self):
        return f"{self.b_seq}-{self.y_seq}"

    @property
    def fasta_name(self):
        return f"{HS_PREFIX}{self.b_seq}-{self.y_seq}"

    @classmethod
    def from_name(cls, name: str):
        b_seq, y_seq = cls.parse_name_to_b_and_y_seqs(name=name)
        return cls(
            b_seq=b_seq,
            y_seq=y_seq,
        )

    @staticmethod
    def parse_name_to_b_and_y_seqs(name: str):
        """
        Parse the name of a hybrid peptide to extract the b and y sequences.
        Depends on how .fasta_name works
        """
        if name.startswith(HS_PREFIX):
            name = name[len(HS_PREFIX) :]
        # For back compatibility
        elif name.startswith("hybrid_"):
            name = name[len("hybrid_") :]
        else:
            raise ValueError(
                f"Invalid hybrid peptide name: {name}. Expected prefix '{HS_PREFIX}' or 'hybrid_'."
            )

        b_seq, y_seq = name.split("-")
        return b_seq, y_seq

    def set_protein_names(self, prot_id_to_name_map: Dict[int, str]):
        self.b_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.b_prot_ids
        )
        self.y_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.y_prot_ids
        )

    def set_fasta_info(
        self,
        prot_id_to_name_map: Optional[Dict[int, str]],
    ):
        self.set_protein_names(prot_id_to_name_map=prot_id_to_name_map)
        self.fasta_description = f"b-prots:{','.join(self.b_prot_names)} y-prots:{','.join(self.y_prot_names)}"

    def mz(self, charge: int):
        return compute_peptide_mz(aa_seq=self.b_seq + self.y_seq, charge=charge)

    def to_dict(self):
        return {
            "b_seq": self.b_seq,
            "y_seq": self.y_seq,
            "b_prot_ids": list(self.b_prot_ids),
            "y_prot_ids": list(self.y_prot_ids),
            "b_prot_names": list(self.b_prot_names) if self.b_prot_names else None,
            "y_prot_names": list(self.y_prot_names) if self.y_prot_names else None,
        }

    @classmethod
    def from_json(cls, json: Path) -> List["HybridPeptide"]:
        hybrid_seq_to_dicts = load_json(in_path=json)
        hybrids = []
        for _, hybrid_dicts in hybrid_seq_to_dicts.items():
            for hybrid_dict in hybrid_dicts:
                hybrid_dict["sample"] = json.parent.stem
                hybrid_dict["scan"] = int(json.stem)
                hybrids.append(cls(**hybrid_dict))
        return hybrids


@dataclass
class SeqWithMass(SqlTableRow):
    seq: str
    mz: float

    @classmethod
    def from_seq(cls, seq: str, charge: int):
        mz = compute_peptide_mz(aa_seq=seq, charge=charge)
        return cls(seq=seq, mz=mz)


def hypedsearch_output_stem(mzml_stem: str):
    return f"{HS_PREFIX}{mzml_stem}"


@log_time(level=logging.DEBUG)
def remove_native_hybrids(
    proteins: List[Peptide],
    hybrids: List[HybridPeptide],
    min_k: int,
    max_k: int,
) -> List[HybridPeptide]:
    kmer_to_prot_id_map = get_uniq_kmer_to_protein_map(
        proteins=proteins, min_k=min_k, max_k=max_k
    )
    non_native_hybrids = list(
        filter(lambda hybrid: hybrid.seq not in kmer_to_prot_id_map, hybrids)
    )
    return non_native_hybrids


def parse_hs_result_file_name(file_name: str):
    hs_stem, scan, ext = file_name.split(".")
    scan_num = int(scan.split("-")[0])
    mzml_stem = hs_stem[len(HS_PREFIX) :]
    return mzml_stem, scan_num


def hybrid_fasta_name(hybrid_seq: str) -> str:
    return f"{HS_PREFIX}{hybrid_seq}"


def create_hybrids_fasta(
    hybrid_seqs: List[str],
    new_fasta_path: Path,
    old_fasta: Optional[Path] = None,
    protein_names: Optional[Union[List[str], Path]] = None,
) -> List[Peptide]:
    """
    Writes a list of hybrids (hybrids) to a FASTA file (new_fasta_path).
    If old_fasta is provided, it will also include the proteins from the old FASTA.
    If protein_names is provided, it will only include those proteins from the old FASTA.
    """
    prots = []
    if old_fasta is not None:
        if protein_names is not None:
            # Get specific proteins from FASTA by name
            prots = get_proteins_by_name(
                protein_names=protein_names, fasta_path=old_fasta
            )
        else:
            # Get all the proteins in the FASTA
            prots = Peptide.from_fasta(fasta_path=old_fasta)

    for _, hybrid_seq in enumerate(hybrid_seqs):
        new_peptide = Peptide(
            seq=hybrid_seq,
            name=hybrid_fasta_name(hybrid_seq=hybrid_seq),
        )

        prots.append(new_peptide)

    write_fasta(peptides=prots, output_path=new_fasta_path)
    return prots


@log_time(level=logging.DEBUG)
def postprocess_hybrids(
    hybrids: List[HybridPeptide],
    fasta_path: Path,
    prot_id_to_name_map: Dict[int, str],
    remove_native: bool = True,
) -> Dict[str, List[HybridPeptide]]:
    """
    Remove hybrids that are native sequences (i.e., they are in the FASTA)
    and group hybrids by their sequence (e.g., group A-BC with AB-C)
    """
    logger.debug(
        f"Postprocessing hybrids... Initially there are {len(hybrids)} hybrids"
    )
    fasta = Fasta(fasta_path=fasta_path)

    seq_to_hybrid_peptides = defaultdict(list)
    for hybrid in hybrids:
        hybrid.set_protein_names(prot_id_to_name_map=prot_id_to_name_map)
        seq_to_hybrid_peptides[hybrid.seq].append(hybrid)

    # Remove hybrids that are native sequences
    if remove_native:
        seq_to_hybrid_peptides = {
            seq: hybrids
            for seq, hybrids in seq_to_hybrid_peptides.items()
            if not fasta.contains_seq(query_seq=seq)
        }
    logger.debug(
        f"Finished postprocessing hybrids. There are {len(seq_to_hybrid_peptides)} hybrid sequences remaining"
    )
    return seq_to_hybrid_peptides
