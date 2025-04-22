# makes all fields keyword-only so ordering nonsense (like "TypeError: non-default
# argument 'scan' follows default argument") doesn't happen
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Dict, List, Optional, Union

import pandas as pd

from src.comet_utils import CometPSM
from src.constants import GIT_REPO_DIR
from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.utils import flatten_list_of_lists


@dataclass(kw_only=True)
class TrueHybrid(HybridPeptide):
    sample: ClassVar[str] = "BMEM_AspN_Fxn5"
    psms: List[CometPSM] = field(init=False)
    hybrid_containing_psms: List[CometPSM] = field(init=False)
    in_psms: bool = field(init=False)
    scan: int
    spectrum: Optional[Spectrum] = None
    name_mapping: ClassVar[Dict] = {
        "Scg1": [
            "sp|P16014|SCG1_MOUSE",
            "sp|Q03517|SCG2_MOUSE",
            "sp|P47867|SCG3_MOUSE",
        ],
        "Ins1": ["sp|P01325|INS1_MOUSE"],
        "Ins2": ["sp|P01326|INS2_MOUSE"],
        "Scg": ["sp|Q03517|SCG1_MOUSE"],
        "IAPP": ["sp|P12968|IAPP_MOUSE", "sp|P12968B|IAPP_MOUSE"],
    }

    def __repr__(self):
        # Only include fields that are not None or unset
        set_fields = {
            key: value for key, value in vars(self).items() if value is not None
        }
        # Creating a string representation with just the set attributes
        return f"TrueHybrid({', '.join(f'{k}={v}' for k, v in set_fields.items())})"

    def get_spectrum(self) -> Spectrum:
        self.spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=self.sample, scan_num=self.scan
        )
        return self.spectrum

    @classmethod
    def parse_protein_column(cls, prot: str) -> List[str]:
        # There can be values like "Scg1" or "Ins1/Ins2"

        if "/" in prot:
            prots = prot.split("/")
        else:
            prots = [prot]
        return flatten_list_of_lists([cls.name_mapping[prot] for prot in prots])

    @classmethod
    def from_df_row(cls, row: pd.Series):
        b_seq, y_seq = row["Sequence"].split("-")
        return cls(
            scan=row["Scan Number"],
            b_seq=b_seq,
            y_seq=y_seq,
            b_prot_names=set(cls.parse_protein_column(prot=row["Left Protein"])),
            y_prot_names=set(cls.parse_protein_column(prot=row["Right protein"])),
        )

    def compare_to_comet(self, comet_txt: Path) -> List[CometPSM]:
        self.psms = CometPSM.from_txt(file_path=comet_txt, sample=self.sample)
        self.in_psms = seq_in_psms(seq=self.seq, psms=self.psms)
        self.hybrid_containing_psms = get_seq_containing_psms(
            seq=self.seq, psms=self.psms
        )

    def set_hybrid_in_psms(self) -> bool:
        psms_containing_seq = list(filter(lambda psm: self.seq in psm.seq, self.psms))
        if len(psms_containing_seq) > 0:
            self.hybrid_in_psms = True
        else:
            self.hybrid_in_psms = False

        return self.hybrid_in_psms

    @property
    def highest_scoring_psm(self):
        return get_highest_scoring_psm(psms=self.psms)


@dataclass(kw_only=True)
class AllHybridsRun:
    hybrid: TrueHybrid
    comet_txt: Path

    def compare(self):
        self.hybrid.compare_to_comet(comet_txt=self.comet_txt)


def get_seq_containing_psms(seq: str, psms: List[CometPSM]) -> List[CometPSM]:
    psms_containing_seq = list(filter(lambda psm: seq in psm.seq, psms))
    return psms_containing_seq


def seq_in_psms(seq: str, psms: List[CometPSM]):
    psms_containing_seq = get_seq_containing_psms(seq=seq, psms=psms)

    if len(psms_containing_seq) > 0:
        return True
    else:
        return False


def get_highest_scoring_psm(psms: List[CometPSM]) -> CometPSM:
    return max(psms, key=lambda psm: psm.xcorr)


def get_true_hybrids_from_thomas_file(
    file_path: Union[str, Path] = GIT_REPO_DIR / "data/found_hybrids.csv",
):
    # Read Thomas's file
    true_hybrids_df = pd.read_csv(file_path)

    # Get the Thomas names I need mappings for
    print(
        f"Protein names in true hybrids file: {set(true_hybrids_df['Right protein'].unique()).union(true_hybrids_df['Left Protein'].unique())}"
    )

    # Create objects better represent the true hybrids
    true_hybrids = [
        TrueHybrid.from_df_row(row=row) for _, row in true_hybrids_df.iterrows()
    ]
    return true_hybrids
