from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import pandas as pd

from src.comet_utils import read_comet_txt_to_df
from src.constants import HS_OUTPUT_PREFIX, HybridFormingMethods
from src.peptides_and_ions import (
    Peptide,
    compute_peptide_mz,
    get_uniq_kmer_to_protein_map,
    write_fasta,
)
from src.sql_database import SqlTableRow
from src.utils import log_time


@dataclass
class HybridPeptide:
    b_seq: str
    y_seq: str
    b_prot_ids: Set[int] = field(default_factory=set)
    b_prot_names: Set[str] = field(default_factory=set)
    y_prot_ids: Set[int] = field(default_factory=set)
    y_prot_names: Set[str] = field(default_factory=set)

    @property
    def seq(self):
        return self.b_seq + self.y_seq

    @property
    def seq_with_hyphen(self):
        return f"{self.b_seq}-{self.y_seq}"

    @property
    def fasta_name(self):
        return f"hybrid_{self.b_seq}-{self.y_seq}"

    def set_protein_names(self, prot_id_to_name_map: Dict[int, str]):
        self.b_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.b_prot_ids
        )
        self.y_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.y_prot_ids
        )

    def set_fasta_info(
        self,
        prot_id_to_name_map: Dict[int, str],
    ):
        self.set_protein_names(prot_id_to_name_map=prot_id_to_name_map)
        self.fasta_description = f"b-prots:{','.join(self.b_prot_names)} y-prots:{','.join(self.y_prot_names)}"

    def mz(self, charge: int):
        return compute_peptide_mz(aa_seq=self.b_seq + self.y_seq, charge=charge)


@dataclass
class HSOutput:
    hybrids: List[HybridPeptide]
    comet_txt: Optional[Path]
    run_time: Optional[float] = None


def hypedsearch_output_stem(mzml_stem: str):
    return f"{HS_OUTPUT_PREFIX}{mzml_stem}"


def determe_hybrid_forming_method(
    protein_names: Optional[Path],
) -> HybridFormingMethods:
    if protein_names is not None:
        return HybridFormingMethods.all
    else:
        raise RuntimeError("Unrecognized hybrid forming method!")


@log_time
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
    mzml_stem = hs_stem[len(HS_OUTPUT_PREFIX) :]
    return mzml_stem, scan_num


def get_all_hybrid_run_psms_drom_dir(folder: Path):
    run_result_txts = list(folder.glob(f"{HS_OUTPUT_PREFIX}*"))
    comet_psms = []
    for comet_txt in run_result_txts:
        # sample, _ = parse_hs_result_file_name(file_name=comet_txt.name)
        comet_psms.append(read_comet_txt_to_df(txt_path=comet_txt))
    comet_psms = pd.concat(comet_psms, ignore_index=True)
    return comet_psms


def create_hybrids_fasta(
    hybrids: List[HybridPeptide],
    fasta_path: Path,
    other_prots: Optional[Union[List[Peptide], Path]] = None,
) -> List[Peptide]:
    """
    Writes a list of hybrids (hybrids) write to a FASTA file (fasta_path).
    You can optionally include other proteins in the FASTA file with the hybrids.
    You can include other proteins via other_prots in three different ways: if other_prots
     is a
        1) list - you passed a list of Peptide objects
        2) None - there are no other proteins/peptides you want in the FASTA file; just the hybrids
        3) Path - you passed a FASTA file; all proteins in the FASTA file will be included
         with the hybrids
    """
    prots = []
    if other_prots is not None:
        if isinstance(other_prots, list):
            prots = other_prots
        elif isinstance(other_prots, Path):
            prots = Peptide.from_fasta(fasta_path=other_prots)
    for idx, hybrid in enumerate(hybrids):
        new_peptide = Peptide(
            seq=hybrid.seq, name=hybrid.fasta_name, desc=hybrid.fasta_description
        )

        prots.append(new_peptide)

    write_fasta(peptides=prots, output_path=fasta_path)
    return prots


@dataclass
class SeqWithMass(SqlTableRow):
    seq: str
    mz: float

    @classmethod
    def from_seq(cls, seq: str, charge: int):
        mz = compute_peptide_mz(aa_seq=seq, charge=charge)
        return cls(seq=seq, mz=mz)
