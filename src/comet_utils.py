import logging
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from functools import cached_property
from itertools import groupby
from pathlib import Path
from typing import (
    Callable,
    Counter,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    Self,
    Set,
    Tuple,
    Union,
)

import pandas as pd
from pydantic import BaseModel
from scipy.stats import percentileofscore

from src.constants import (
    COMET,
    CRUX,
    DELTA_CN,
    EVAL,
    HS_PREFIX,
    IONS_MATCHED,
    IONS_TOTAL,
    NUM,
    PLAIN_PEPTIDE,
    PROTEIN,
    Q_VALUE,
    SAMPLE,
    SCAN,
    XCORR,
)
from src.kmer_database import KmerToProteinsMap
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Peptide
from src.utils import (
    flatten_list_of_lists,
    get_arg_fcn_of_objects,
    get_fcn_of_objects,
    load_json,
)

logger = logging.getLogger(__name__)


@dataclass
class CometTxt:
    path: Union[str, Path]
    sample: str = field(init=False)
    file_type: Literal["comet", "crux"] = field(init=False)
    # hybrid_run: Optional[bool]

    def __post_init__(self):
        # Set file_type
        with open(self.path, "r") as f:
            first_line = f.readline().strip()
            if ("comet" in first_line) or ("Comet" in first_line):
                self.file_type = COMET
            else:
                self.file_type = CRUX

        # Set sample
        self.sample = self.path.stem.split(".")[0]

    def read_psms(self, as_df: bool = False) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads the Comet output file to a list of dataclasses or a dataframe
        """
        if self.path is not None:
            return CometPSM.from_txt(txt=self.path, as_df=as_df)
        else:
            raise ValueError("No Comet output file found!")

    def get_header(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[1]

    def get_first_psm_line(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[2]

    def get_top_psms(self) -> List["CometPSM"]:
        psms = CometPSM.from_txt(txt=self.path)
        return [psm for psm in psms if psm.num == 1]


class HybridPeptide(BaseModel):
    left_seq: str
    right_seq: str
    left_proteins: List[str]
    right_proteins: List[str]

    @property
    def seq(self) -> str:
        return self.left_seq + self.right_seq

    def get_left_protein_max_percentile(self, protein_counts: Counter) -> float:
        return max(
            [
                percentileofscore(list(protein_counts.values()), protein_counts[prot])
                for prot in self.left_proteins
            ],
        )

    def get_right_protein_max_percentile(self, protein_counts: Counter) -> float:
        return max(
            [
                percentileofscore(list(protein_counts.values()), protein_counts[prot])
                for prot in self.right_proteins
            ],
        )

    def get_left_seq_kmer_percentile(self, kmer_counts: Dict[int, Counter]) -> float:
        return percentileofscore(
            list(kmer_counts[len(self.left_seq)].values()),
            kmer_counts[len(self.left_seq)][self.left_seq],
        )

    def get_right_seq_kmer_percentile(self, kmer_counts: Dict[int, Counter]) -> float:
        return percentileofscore(
            list(kmer_counts[len(self.right_seq)].values()),
            kmer_counts[len(self.right_seq)][self.right_seq],
        )


def read_comet_psms_from_dir(
    dir_path: Union[str, Path], glob_pattern: Optional[str] = None
):
    if glob_pattern is None:
        glob_pattern = "*.txt"
    all_psms = []
    for txt_file in Path(dir_path).glob(glob_pattern):
        psms = CometPSM.from_txt(txt=txt_file)
        all_psms.extend(psms)
    return all_psms


@dataclass
class CometPSM:
    """Class for rows of Comet output"""

    sample: str
    num: int
    scan: int
    seq: str
    ions_matched: int
    ions_total: int
    proteins: List[str]
    # protein_count: int
    xcorr: float
    eval: float
    delta_cn: float
    q_value: Optional[float]

    @classmethod
    def from_txt(
        cls,
        txt: str,
        as_df: bool = False,
        sample: str = "",
    ) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads Comet results .txt file to a list of dataclasses or a dataframe
        """
        # Check whether TXT is from a direct Comet run or a Comet run via crux
        comet_txt = CometTxt(path=Path(txt))

        # Set sample if not provided
        if len(sample) == 0:
            sample = comet_txt.sample
        if comet_txt.file_type == CRUX:
            df = pd.read_csv(comet_txt.path, sep="\t")
            df[SAMPLE] = sample
            if "file" in df.columns:
                # If the 'file' column exists, it means it's the output of `crux assign-confidence`
                # in which case we need to set the sample differently
                df[SAMPLE] = df["file"].apply(
                    lambda file_path: Path(file_path).stem.split(".")[0]
                )
            df.rename(
                columns={
                    "b/y ions matched": IONS_MATCHED,
                    "b/y ions total": IONS_TOTAL,
                    "xcorr score": XCORR,
                    "xcorr rank": NUM,
                    "protein id": PROTEIN,
                    "sequence": PLAIN_PEPTIDE,
                    "tdc q-value": Q_VALUE,
                },
                inplace=True,
            )
        elif comet_txt.file_type == COMET:
            df = pd.read_csv(comet_txt.path, sep="\t", header=1)
            df[SAMPLE] = sample

        if as_df:
            return df
        else:
            return [
                cls(
                    sample=row[SAMPLE],
                    scan=row[SCAN],
                    num=row[NUM],
                    ions_matched=row[IONS_MATCHED],
                    ions_total=row[IONS_TOTAL],
                    # protein_count=row[PROTEIN_COUNT],
                    proteins=row[PROTEIN].split(","),
                    seq=row[PLAIN_PEPTIDE],
                    xcorr=row[XCORR],
                    eval=row[EVAL],
                    delta_cn=row[DELTA_CN],
                    q_value=row.get(Q_VALUE, None),  # Handle optional q-value
                )
                for _, row in df.iterrows()
            ]

    @property
    def seq_with_hyphen(self):
        if self.is_hybrid:
            hybrid_peptides = self.get_hybrid_peptides()
            return [f"{pep.b_seq}-{pep.y_seq}" for pep in hybrid_peptides]
        else:
            return self.seq

    @staticmethod
    def check_if_hybrid_prot(prot: str):
        if prot.startswith(HS_PREFIX) or prot.startswith("hybrid_"):
            return True
        else:
            return False

    @property
    def is_hybrid(self):
        """
        A Comet PSM is a hybrid if the only proteins it appears in are hybrid proteins.
        If a PSM is in both a hybrid protein and a native protein, that means that the "hybrid"
        is a native sequence.
        """
        if all(self.check_if_hybrid_prot(prot=prot) for prot in self.proteins):
            return True
        else:
            return False

    def remake(self) -> "CometPSM":
        """
        Remakes the CometPSM object from the original PSM
        """
        return CometPSM(
            sample=self.sample,
            num=self.num,
            scan=self.scan,
            seq=self.seq,
            ions_matched=self.ions_matched,
            proteins=self.proteins,
            protein_count=self.protein_count,
            xcorr=self.xcorr,
            eval=self.eval,
            delta_cn=self.delta_cn,
        )

    def to_dict(self):
        return asdict(self)

    def get_possible_hybrid_peptides(
        self, kmer_to_proteins_map: Dict[str, List[str]], min_side_len: int
    ) -> List[HybridPeptide]:
        return get_possible_hybrid_peptides_for_seq(
            seq=self.seq,
            kmer_to_proteins_map=kmer_to_proteins_map,
            min_side_len=min_side_len,
        )

    @property
    def mzml_name(self) -> str:
        return Mzml.get_mzml_name(mzml=self.sample)


def get_possible_hybrid_peptides_for_seq(
    seq: str, kmer_to_proteins_map: Dict[str, List[str]], min_side_len: int
) -> List[HybridPeptide]:
    possible_hybrids = []
    for breakpoint in range(min_side_len, len(seq) - min_side_len + 1):
        left = seq[:breakpoint]
        right = seq[breakpoint:]
        if (left in kmer_to_proteins_map) and (right in kmer_to_proteins_map):
            possible_hybrids.append(
                HybridPeptide(
                    left_seq=left,
                    right_seq=right,
                    left_proteins=kmer_to_proteins_map[left],
                    right_proteins=kmer_to_proteins_map[right],
                )
            )
    return possible_hybrids
