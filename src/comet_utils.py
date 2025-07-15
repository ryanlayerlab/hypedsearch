import logging
import os
import platform
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from functools import cached_property
from itertools import groupby
from pathlib import Path
from typing import Annotated, Callable, List, Literal, Optional, Tuple, Union

import click
import pandas as pd
from pydantic import BaseModel, BeforeValidator, validator

from src.constants import (
    COMET,
    COMET_DIR,
    COMET_PARAMS,
    COMET_RUN_1_DIR,
    COMET_RUN_2_DIR,
    CRUX,
    DEFAULT_COMET_PARAMS_FILE,
    DEFAULT_COMET_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_COMET_SCAN_RANGE,
    DEFAULT_NUM_PSMS,
    DEFAULT_PPM_TOLERANCE,
    DELTA_CN,
    EVAL,
    HS_PREFIX,
    IONS_MATCHED,
    IONS_TOTAL,
    MOUSE_PROTEOME,
    NUM,
    PLAIN_PEPTIDE,
    PROTEIN,
    PROTEIN_COUNT,
    Q_VALUE,
    SAMPLE,
    SCAN,
    THOMAS_SAMPLES,
    XCORR,
)
from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.utils import (
    ExistingPath,
    PathType,
    flatten_list_of_lists,
    get_arg_fcn_of_objects,
    get_default_comet_executable_path,
    get_fcn_of_objects,
    get_os,
    log_params,
    log_time,
    make_directory,
    remove_gene_name,
    setup_logger,
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
            return CometPSM.from_txt(txt_path=self.path, as_df=as_df)
        else:
            raise ValueError("No Comet output file found!")

    def get_header(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[1]

    def get_first_psm_line(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[2]

    def get_top_psms(self) -> List["CometPSM"]:
        psms = CometPSM.from_txt(txt_path=self.path)
        return [psm for psm in psms if psm.num == 1]


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
        txt_path: str,
        as_df: bool = False,
        sample: str = "",
    ) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads Comet results .txt file to a list of dataclasses or a dataframe
        """
        # Check whether TXT is from a direct Comet run or a Comet run via crux
        comet_txt = CometTxt(path=Path(txt_path))

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

    def get_corresponding_spectrum(self) -> Spectrum:
        """
        Gets the spectrum corresponding to the Comet-returned PSM
        """
        return get_specific_spectrum_by_sample_and_scan_num(
            sample=self.sample, scan_num=self.scan
        )

    # @staticmethod
    # def check_if_hybrid(proteins: List[str]) -> bool:
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


@dataclass
class CometPSMs:
    psms: List[CometPSM]

    def get_psms_for_sample_and_scan(self, sample: str, scan: int):
        return list(
            filter(lambda psm: (psm.sample == sample) and (psm.scan == scan), self.psms)
        )

    def loop_over_psms_by_sample_and_scan(self):
        self.psms = sorted(self.psms, key=lambda x: (x.sample, x.scan))
        for key, group in groupby(self.psms, key=lambda x: (x.sample, x.scan)):
            yield key[0], key[1], list(group)

    def get_fcn_of_psms(self, attr: str, fcn: Callable):
        return get_fcn_of_objects(objs=self.psms, attr=attr, fcn=fcn)

    def get_arg_fcn_of_psms(self, attr: str, fcn: Callable):
        return get_arg_fcn_of_objects(objs=self.psms, attr=attr, fcn=fcn)


@dataclass
class ScanPSMs:
    native_psms: CometPSMs
    hybrid_psms: CometPSMs

    def __post_init__(self):
        # Handle if psms are passed as List[CometPSM] instead of CometPSMs
        if isinstance(self.native_psms, list):
            self.native_psms = CometPSMs(psms=self.native_psms)
        if isinstance(self.hybrid_psms, list):
            self.hybrid_psms = CometPSMs(psms=self.hybrid_psms)

    @property
    def sample(self):
        return self.native_psms.psms[0].sample

    @property
    def scan(self):
        return self.native_psms.psms[0].scan


def get_comet_protein_counts(
    comet_rows: Optional[List[CometPSM]] = None,
    shorten_names: bool = True,
) -> Counter:
    # if comet_rows is None:
    #     comet_rows = load_comet_data(as_df=False)

    all_prots = flatten_list_of_lists([row.proteins for row in comet_rows])

    if shorten_names is True:
        all_prots = [remove_gene_name(protein_name=prot) for prot in all_prots]

    return Counter(all_prots)
