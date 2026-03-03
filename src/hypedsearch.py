import logging
import multiprocessing as mp
import os
import re
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from functools import cached_property, partial
from pathlib import Path
from typing import Callable, Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import numpy as np
import pandas as pd
from fm_index import MultiFMIndex
from matplotlib.axes import Axes
from pydantic import BaseModel, ConfigDict, field_validator, model_validator
from scipy.interpolate import PchipInterpolator
from typing_extensions import Self

from src.constants import (
    ASSIGN_CONFIDENCE,
    DECOY,
    DEFAULT_MAX_ALLOWED_ION_CHARGE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MAX_PRECURSOR_CHARGE,
    DEFAULT_MIN_CLUSTER_LENGTH,
    DEFAULT_MIN_CLUSTER_SUPPORT,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_Q_THRESHOLD,
    HUMAN_PROTEOME,
    MAC_CRUX_EXECUTABLE,
    NAT_DECOY,
    NAT_TARGET,
    Q_VAL,
    SHARED_PARAMS,
    TARGET,
    TRUE_HYBRIDS_PATH,
    XCORR,
)
from src.crux import (
    CometOutputs,
    CometRun,
    Crux,
    get_expected_comet_outputs_for_mzml_to_scans,
)
from src.hybrids_via_clusters import HybridPeptide, form_spectrum_hybrids_via_clustering
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum, create_sample_scan_to_spectrum_map
from src.peptides_and_ions import Fasta, Fasta2MFMIndex, Peptide, get_proteins_by_name
from src.plot_utils import (
    fig_setup,
    finalize,
    save_fig,
    score_histogram,
    set_title_axes_labels,
)
from src.psm import CometPSM, ProteinAbundance
from src.utils import (
    PathType,
    flatten_list_of_lists,
    from_pickle,
    load_json,
    log_params,
    mass_difference_in_ppm,
    path_aware_dict_factory,
    save_dict,
    save_pydantic_objects_to_json,
    setup_logger,
    to_pickle,
)

logger = logging.getLogger(__name__)


class TrueHybrid(HybridPeptide):
    spectra: List[Spectrum] = field(default_factory=list)
    rt: Optional[float] = None  # retention time
    id: Optional[Union[str, int]] = None
    precursor_mz: Optional[float] = None

    @classmethod
    def from_excel(
        cls,
        path: Union[str, Path],
        hyphen_seq_colm: str = "hyphen_seq",
        min_side_len: Optional[int] = None,
    ) -> List["TrueHybrid"]:
        df = pd.read_excel(path)

        hybrids = []
        for row_idx, row in df.iterrows():
            try:
                left_seq, right_seq = row[hyphen_seq_colm].split("-")
            except:
                logger.info(
                    f"Skipping sequence {row[hyphen_seq_colm]} because it couldn't be split by `-`"
                )
                continue
            if min_side_len is not None:
                if (len(left_seq) < min_side_len) or (len(right_seq) < min_side_len):
                    logger.info(
                        f"Hybrid {left_seq}-{right_seq} because one side is shorter than min_side_len={min_side_len}"
                    )
                    continue
            hybrids.append(
                cls(
                    id=row.id,
                    left_seq=left_seq,
                    right_seq=right_seq,
                    rt=row.rt if isinstance(row.rt, (int, float)) else None,
                    precursor_mz=row.mz if isinstance(row.mz, (int, float)) else None,
                )
            )
        return hybrids

    @staticmethod
    def save(
        true_hybrids: List["TrueHybrid"],
        out: Union[Path, str] = TRUE_HYBRIDS_PATH,
    ):
        save_pydantic_objects_to_json(objects=true_hybrids, path=out)

    @classmethod
    def load(cls, path: Union[Path, str] = TRUE_HYBRIDS_PATH) -> List["TrueHybrid"]:
        data = load_json(path=path)
        return [cls.model_validate(item) for item in data]

    @staticmethod
    def get_spectra_for_true_hybrids(
        true_hybrids: List["TrueHybrid"],
        mzmls: List[Path],
        precursor_mz_ppm_tol: float,
        retention_time_tol: float,
    ) -> None:
        sample_scan_to_spectrum_map = create_sample_scan_to_spectrum_map(mzmls=mzmls)
        for sample_scan, spectrum in sample_scan_to_spectrum_map.items():
            for hybrid in true_hybrids:
                try:
                    ppm_diff = mass_difference_in_ppm(
                        mass1=spectrum.precursor_mz, mass2=hybrid.precursor_mz
                    )
                    rt_diff = abs(spectrum.retention_time - hybrid.rt)
                except:
                    continue
                if ppm_diff <= precursor_mz_ppm_tol and rt_diff <= retention_time_tol:
                    hybrid.spectra.append(spectrum)

    @staticmethod
    def add_proteins(
        true_hybrids: List["TrueHybrid"],
        fasta: Union[Fasta, Path, str],
    ) -> List["TrueHybrid"]:
        if isinstance(fasta, (str, Path)):
            fasta = Fasta(path=fasta)
        seq_to_prots = fasta.proteins_that_contain_seqs(
            seqs=flatten_list_of_lists(
                [[hy.left_seq, hy.right_seq] for hy in true_hybrids]
            )
        )
        for hy in true_hybrids:
            hy.left_proteins = seq_to_prots[hy.left_seq]
            hy.right_proteins = seq_to_prots[hy.right_seq]
        return true_hybrids


@dataclass
class HypedsearchOutputs:
    comet_regex = r"^(?P<mzml>(.+?))\.comet\.(?P<scan>\d+)-(?P=scan)\.(?P<psm_type>target|decoy)\.txt$"
    assign_conf_regex = r"^(?P<name>(.+?))-assign-confidence.txt"

    @staticmethod
    def get_combined_comet_output_name(
        mzml_name: str, psm_type: Literal[DECOY, TARGET]
    ) -> str:
        return f"{mzml_name}.comet.{psm_type}.txt"

    @staticmethod
    def get_expected_output_txt_name(
        name: str,
        psm_type: Literal[TARGET, DECOY, ASSIGN_CONFIDENCE],
        scan: Optional[int] = None,
    ):
        if psm_type == ASSIGN_CONFIDENCE:
            return f"{name}-{ASSIGN_CONFIDENCE}.txt"
        else:
            return f"{name}.comet.{scan}-{scan}.{psm_type}.txt"

    @classmethod
    def get_comet_outputs_in_folder(
        cls,
        folder: Union[str, Path],
    ) -> Dict[str, Dict[str, Path]]:
        """
        Given a folder, find all Comet output txt files and group them by MZML name and PSM type (=target or decoy). Return a dictionary of the form
        """
        # Get Comet outputs for each MZML
        logger.info("Grouping Comet outputs by MZML...")
        mzml_to_psm_type_to_txts = defaultdict(lambda: {TARGET: [], DECOY: []})
        for comet_txt in Path(folder).glob("*.txt"):
            match = re.match(cls.comet_regex, comet_txt.name)
            if match is None:
                logger.info(
                    f"Skipping file {comet_txt} because it doesn't match expected pattern"
                )
                continue

            mzml_name = match.groupdict()["mzml"]
            psm_type = match.groupdict()["psm_type"]
            mzml_to_psm_type_to_txts[mzml_name][psm_type].append(comet_txt)
        return dict(mzml_to_psm_type_to_txts)

    @classmethod
    def parse_combined_comet_output_name(cls, filename: str):
        match = re.match(cls.comet_regex, filename)
        if match is not None:
            mzml_name = match.groupdict()["mzml"]
            psm_type = match.groupdict()["psm_type"]
            return mzml_name, psm_type

        match = re.match(cls.assign_conf_regex, filename)
        if match is not None:
            mzml_name = match.groupdict()["name"]
            psm_type = ASSIGN_CONFIDENCE
            return mzml_name, psm_type

        else:
            raise ValueError(
                f"Filename {filename} does not match expected combined comet output patterns: {cls.assign_conf_regex} nor {cls.comet_regex}"
            )

    @classmethod
    def combine_comet_results_by_mzml_and_psm_type(
        cls, folder: Union[Path, str], out_dir: Optional[Union[Path, str]] = None
    ) -> Dict[str, Dict[str, Path]]:
        if out_dir is None:
            out_dir = Path(folder).parent
            logger.info(f"Out directory not set so setting to {out_dir}")
        # Get Comet outputs for each MZML
        mzml_to_psm_type_to_txts = cls.get_comet_outputs_in_folder(folder=folder)
        for mzml_name, psm_type_to_txts in mzml_to_psm_type_to_txts.items():
            for psm_type, txts in psm_type_to_txts.items():
                if len(txts) == 0:
                    continue
                logger.info(f"Combining Comet {psm_type} outputs for {mzml_name}...")
                _ = Crux.combine_crux_comet_files(
                    files=txts,
                    out_path=Path(out_dir)
                    / cls.get_combined_comet_output_name(
                        mzml_name=mzml_name, psm_type=psm_type
                    ),
                )

    @classmethod
    def get_matching_output_txt(out_dir: Path, mzml_name: str, psm_type: str) -> Path:
        path
        for path in out_dir.iterdir():
            if path.is_file():
                try:
                    found_mzml_name, found_psm_type = parse_combined_comet_output_name(
                        filename=path.name
                    )
                    if found_mzml_name == mzml_name and found_psm_type == psm_type:
                        return path
                except:
                    continue
        raise RuntimeError(
            f"No matching output txt found for MZML={mzml_name}, PSM_TYPE={psm_type} in directory {out_dir}"
        )


class HybridPSMScorer(BaseModel):
    comet_params: Path
    hybrid_decoy_competition: bool = False
    fasta: Optional[Path] = None

    @model_validator(mode="after")
    def post_init(self) -> Self:
        if self.hybrid_decoy_competition and self.fasta is None:
            raise ValueError(
                "A FASTA must be provided if hybrid_decoy_competition is True"
            )
        return self

    @property
    def decoy_search(self) -> Literal[0, 2]:
        if self.hybrid_decoy_competition:
            return 2
        else:
            return 0

    def score_hybrids(
        self,
        seq_to_hybrids: Dict[str, List[HybridPeptide]],
        spectrum: Spectrum,
        out_dir: Path,
    ) -> CometOutputs:
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create FASTA containing hybrids and run Comet
            tmp_path = Path(tmp_dir)
            mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
            hybrids_fasta = tmp_path / f"{mzml_name}.{spectrum.scan}.fasta"
            if self.hybrid_decoy_competition:
                create_hybrids_fasta(
                    seq_to_hybrids=seq_to_hybrids,
                    output_fasta_path=hybrids_fasta,
                    fasta_to_include=self.fasta,
                )
            else:
                create_hybrids_fasta(
                    seq_to_hybrids=seq_to_hybrids,
                    output_fasta_path=hybrids_fasta,
                )
            # Run Comet
            outputs = Crux().run_comet(
                mzml=spectrum.mzml,
                fasta=hybrids_fasta,
                crux_comet_params=self.comet_params,
                decoy_search=self.decoy_search,
                out_dir=out_dir,
                file_root=mzml_name,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
                num_threads=1,
            )
        return outputs

    def expected_hybrid_target_outputs(
        self, mzml_to_scans: Dict[str, Set[int]], out_dir: Union[Path, str]
    ):
        return get_expected_comet_outputs_for_mzml_to_scans(
            mzml_to_scans=mzml_to_scans,
            out_dir=out_dir,
            decoy_search=self.decoy_search,
            psm_type=TARGET,
        )


class SpectrumSelector(BaseModel):
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE

    def get_selected_spectra(self, spectra: List[Spectrum]) -> List[Spectrum]:
        return [
            spectrum
            for spectrum in spectra
            if spectrum.precursor_charge <= self.max_precursor_charge
        ]

    def get_scan_numbers_of_selected_spectra_from_mzml(self, mzml: Path) -> Set[int]:
        logger.info(f"Getting spectra from {mzml.name}")
        spectra = self.get_selected_spectra(
            spectra=Spectrum.parse_ms2_from_mzml(mzml=mzml)
        )
        logger.info("Done selecting spectra")
        return {spectrum.scan for spectrum in spectra}


class SpectrumPreprocessor(BaseModel):
    num_peaks: int = 0

    def preprocess_spectrum(self, spectrum: Spectrum) -> Spectrum:
        if self.num_peaks > 0:
            logger.info(f"Filtering to top {self.num_peaks} peaks...")
            spectrum.filter_to_top_n_peaks(n=self.num_peaks)
        return spectrum


def get_seq_to_hybrids_map(
    seqs: Union[Set[str], List[str]],
    db_path: Path,
    min_side_len: int = DEFAULT_MIN_SIDE_LEN,
    remove_carbamidomethylation: bool = True,
) -> Dict[str, List[HybridPeptide]]:
    kmer_to_proteins_map = KmerDatabase(
        db_path=db_path
    ).kmer_to_proteins_map.kmer_to_protein_map
    seq_to_hybrids = {}
    for seq in set(seqs):
        hybrids = find_possible_hybrids_for_seq(
            seq=seq,
            kmer_to_proteins_map=kmer_to_proteins_map,
            min_side_len=min_side_len,
        )
        if remove_carbamidomethylation:
            hybrids = [hy for hy in hybrids if not hy.evidence_of_carbamidomethylation]
        if len(hybrids) > 0:
            seq_to_hybrids[seq] = hybrids
    return seq_to_hybrids


def find_possible_hybrids_for_seq(
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


@dataclass
class HybridRunParams:
    kmer_db: Path
    fasta: Path
    crux_comet_params: Path
    fasta_fm_index: MultiFMIndex
    out_dir: Path
    precursor_mz_ppm_tol: float
    peak_to_ion_ppm_tol: float
    min_hybrid_side_len: int
    min_cluster_support: int
    max_allowed_ion_charge: int
    remove_carbamidomethylated_hybrids: bool
    hybrid_decoy_competition: bool


@dataclass
class HypedsearchRunConfig:
    name: str
    mzml_to_scans: Dict[Path, Union[Set[int], Literal["all"]]]
    parent_output_dir: Path
    crux_comet_params: Path
    fasta_fm_index: Path
    kmer_db_path: Optional[Path] = None
    hybrid_decoy_competition: bool = False
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE
    fasta: Path = HUMAN_PROTEOME
    num_peaks: int = 0
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    min_hybrid_side_len: int = DEFAULT_MIN_SIDE_LEN
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE
    remove_carbamidomethylation_hybrids: bool = True

    def __post_init__(self):
        # Set any strings that are supposed to be Path objects to Path objects
        self.mzml_to_scans = {
            Path(mzml): scans for mzml, scans in self.mzml_to_scans.items()
        }
        self.parent_output_dir = Path(self.parent_output_dir)
        self.crux_comet_params = Path(self.crux_comet_params)
        self.fasta_fm_index = Path(self.fasta_fm_index)
        self.fasta = Path(self.fasta)

        # Set kmer database path
        if self.kmer_db_path is None:
            self.kmer_db_path = self.parent_output_dir / f"{self.name}.kmers.db"
        else:
            self.kmer_db_path = Path(self.kmer_db_path)

        # Ensure directories exist
        self.native_run_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_run_scan_results_dir.mkdir(parents=True, exist_ok=True)
        return self

    @property
    def kmer_db(self) -> KmerDatabase:
        return KmerDatabase(db_path=self.kmer_db_path)

    @property
    def mzml_names(self) -> List[str]:
        return [Mzml(path=mzml).name for mzml in self.mzml_to_scans.keys()]

    @property
    def native_run_dir(self) -> Path:
        return self.parent_output_dir / "native_run"

    @property
    def hybrid_run_dir(self) -> Path:
        return self.parent_output_dir / "hybrid_run"

    @property
    def hybrid_run_scan_results_dir(self) -> Path:
        return self.hybrid_run_dir / "scan_results"

    @property
    def spectrum_selector(self) -> SpectrumSelector:
        return SpectrumSelector(max_precursor_charge=self.max_precursor_charge)

    @property
    def psm_scorer(self) -> HybridPSMScorer:
        if self.hybrid_decoy_competition:
            return HybridPSMScorer(
                fasta=self.fasta,
                comet_params=self.crux_comet_params,
                hybrid_decoy_competition=self.hybrid_decoy_competition,
            )
        else:
            return HybridPSMScorer(
                comet_params=self.crux_comet_params,
            )

    @cached_property
    def _mzml_to_spectra(self) -> Dict[Path, List[Spectrum]]:
        logger.info(f"Getting spectra for config {self.name}...")
        mzml_to_spectra = {}
        for mzml, scans in self.mzml_to_scans.items():
            all_selected_spectra = self.spectrum_selector.get_selected_spectra(
                spectra=Spectrum.parse_ms2_from_mzml(mzml=mzml)
            )
            if scans == "all":
                mzml_to_spectra[mzml] = all_selected_spectra
            else:
                # Filter to only spectra with scan numbers in scans
                mzml_to_spectra[mzml] = [
                    spectrum
                    for spectrum in all_selected_spectra
                    if spectrum.scan in scans
                ]
        return mzml_to_spectra

    @cached_property
    def spectrum_uid_to_spectrum(self) -> Dict[str, Spectrum]:
        uid_to_spectrum = {}
        for spectra in self._mzml_to_spectra.values():
            uid_to_spectrum.update({spectrum.uid: spectrum for spectrum in spectra})
        logger.info("Done getting spectra")
        return uid_to_spectrum

    @property
    def spectra(self) -> List[Spectrum]:
        return list(self.spectrum_uid_to_spectrum.values())

    @classmethod
    def from_json(cls, path: Union[Path, str]):
        data = load_json(path=path)
        Path(data["parent_output_dir"]).mkdir(exist_ok=True, parents=True)
        return cls(**data)

    @property
    def native_assign_confidence_path(self) -> Path:
        return self.native_run_dir / HypedsearchOutputs.get_expected_output_txt_name(
            name=self.name, psm_type=ASSIGN_CONFIDENCE
        )

    def run_native_assign_confidence(self):
        native_target_txts = []
        for output in self.expected_native_comet_outputs:
            assert (
                output.target.exists()
            ), f"Expected target txt does not exist: {output.target}"
            native_target_txts.append(output.target)
        Crux().run_assign_confidence(
            target_txts=native_target_txts,
            out_path=self.native_assign_confidence_path,
        )

    @cached_property
    def native_assign_confidence_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.native_assign_confidence_path)

    @property
    def results_dir(self) -> Path:
        return self.parent_output_dir / f"{self.name}_results"

    @property
    def native_results_dir(self) -> Path:
        path = self.parent_output_dir / f"native_{self.name}_results"
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def hybrid_results_dir(self) -> Path:
        path = self.parent_output_dir / f"hybrid_{self.name}_results"
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        return path

    @cached_property
    def _mzml_to_scan_nums(self) -> Dict[Path, List[int]]:
        mzml_to_scan_nums = {}
        for mzml, spectra in self._mzml_to_spectra.items():
            mzml_to_scan_nums[mzml] = [spectrum.scan for spectrum in spectra]
        return mzml_to_scan_nums

    @property
    def expected_hybrid_run_scan_target_txts(self) -> Set[Path]:
        txts = self.psm_scorer.expected_hybrid_target_outputs(
            mzml_to_scans=self._mzml_to_scan_nums,
            out_dir=self.hybrid_run_scan_results_dir,
        )
        assert len(txts) == len(set(txts))
        return set(txts)

    @property
    def existing_hybrid_run_scan_target_txts(self) -> Set[Path]:
        txts = []
        for mzml_name in self.mzml_names:
            txts.extend(
                list(
                    self.hybrid_run_scan_results_dir.glob(
                        f"{mzml_name}.comet.*.{TARGET}.txt"
                    )
                )
            )
        assert len(txts) == len(set(txts))
        return set(txts)

    @property
    def missing_hybrid_run_scan_target_txts(self) -> Set[Path]:
        return (
            self.expected_hybrid_run_scan_target_txts
            - self.existing_hybrid_run_scan_target_txts
        )

    def to_dict(self):
        return asdict(self, dict_factory=path_aware_dict_factory)

    def save(
        self,
        path: Union[Path, str],
    ):
        """Save config as a JSON"""
        save_dict(data=self.to_dict(), path=path)

    def run_native_comet(
        self,
        crux_path: Optional[str] = None,
        on_singularity: bool = False,
        dry_run: bool = False,
    ) -> List[CometOutputs]:
        expected_outputs = []
        mzmls = list(self.mzml_to_scans.keys())
        for idx, mzml in enumerate(mzmls):
            comet_run = CometRun(
                fasta=self.fasta,
                mzml=mzml,
                crux_comet_params=self.crux_comet_params,
                decoy_search=2,
                out_dir=self.native_run_dir,
                file_root=Mzml.get_mzml_name(mzml=mzml),
                dry_run=dry_run,
            )
            expected_outputs.append(comet_run.standardized_comet_outputs)
            if not dry_run:
                logger.info(
                    f"Native Comet run on MZML {mzml.name} ({idx+1}/{len(mzmls)})"
                )
                if comet_run.standardized_comet_outputs.target.exists():
                    logger.info(f"Looks like target TXT already exists so skipping...")
                    continue
                comet_run.run_comet_and_keep_only_results(
                    crux_path=crux_path, on_singularity=on_singularity
                )

        return expected_outputs

    @property
    def expected_native_comet_outputs(self):
        return self.run_native_comet(dry_run=True)

    @cached_property
    def native_target_psms(self):
        logger.info("Getting native target PSMs...")
        txts = []
        for output in self.expected_native_comet_outputs:
            assert (
                output.target.exists()
            ), f"Expected target txt does not exist: {output.target}"
            txts.append(output.target)
        psms = CometPSM.from_txts(txts=txts)
        return psms
        # return self.native_assign_confidence_psms

    @cached_property
    def native_decoy_psms(self):
        logger.info("Getting native decoy PSMs...")
        txts = []
        for output in self.expected_native_comet_outputs:
            assert (
                output.decoy.exists()
            ), f"Expected decoy txt does not exist: {output.decoy}"
            txts.append(output.decoy)
        psms = CometPSM.from_txts(txts=txts)
        return psms

    @cached_property
    def hybrid_target_psms(self):
        logger.info("Getting hybrid target PSMs...")
        txts = []
        for txt in self.expected_hybrid_txts:
            assert txt.exists(), f"Expected hybrid target txt does not exist: {txt}"
            txts.append(txt)
        psms = CometPSM.from_txts(txts=txts)
        return psms

    def get_protein_abundance(
        self,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
    ) -> ProteinAbundance:
        return ProteinAbundance.from_comet_psms(
            quality_psms=self.native_assign_confidence_psms, q_threshold=q_threshold
        )

    def create_kmer_db(
        self,
        top_n_prots: Optional[int] = None,
        min_psm_count: Optional[int] = None,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
        min_k: int = DEFAULT_MIN_KMER_LEN,
        max_k: int = DEFAULT_MAX_KMER_LEN,
    ) -> KmerDatabase:
        # Get proteins
        fasta = Fasta(path=self.fasta)
        prot_cnts = self.get_protein_abundance(
            q_threshold=q_threshold,
        ).protein_counts
        if top_n_prots:
            prots = [prot_cnt[0] for prot_cnt in prot_cnts.most_common(n=top_n_prots)]
        elif min_psm_count:
            prots = [prot for prot, cnt in prot_cnts.items() if cnt >= min_psm_count]
        else:
            raise ValueError("Either top_n_prots or min_psm_count must be provided")
        assert len(prots) > 0, "No proteins selected for kmer database!"
        db = KmerDatabase.create_db(
            db_path=self.kmer_db_path,
            proteins=fasta.get_proteins_by_name(names=prots),
            min_k=min_k,
            max_k=max_k,
            overwrite=True,
        )
        return db

    @property
    def expected_hybrid_txts(self):
        return [
            Path(self.hybrid_run_dir)
            / HypedsearchOutputs.get_combined_comet_output_name(
                mzml_name=mzml_name, psm_type=TARGET
            )
            for mzml_name in self.mzml_names
        ]

    def combine_hybrid_run_scan_outputs(self):
        logger.info(f"Combining Comet scan results for config: {self.name}")
        mzml_to_txts = defaultdict(list)
        for txt in self.expected_hybrid_run_scan_target_txts:
            assert txt.exists(), "Expected txt does not exist: {txt}"
            mzml, _, _ = CometOutputs.parse_standardized_comet_txt(comet_txt=txt)
            mzml_to_txts[mzml].append(txt)
        for mzml, txts in mzml_to_txts.items():
            logger.info(f"Combining Comet scan results for MZML {mzml}...")
            _ = Crux.combine_crux_comet_files(
                files=txts,
                out_path=Path(self.hybrid_run_dir)
                / HypedsearchOutputs.get_combined_comet_output_name(
                    mzml_name=mzml, psm_type=TARGET
                ),
            )

    def analyze_native_results(
        self, q_threshold: float = DEFAULT_Q_THRESHOLD, top_n_prots: int = 100
    ):
        self.create_native_xcorr_plot(save=True)
        prot_ab = self.get_protein_abundance(q_threshold=q_threshold)
        prot_ab.plot_sorted_prot_cnts(top_n_prots=top_n_prots)
        prot_ab.to_json(
            path=self.results_dir / f"protein_abundance_q{q_threshold}.json"
        )

    @property
    def hybrid_run_params(self) -> HybridRunParams:
        return HybridRunParams(
            kmer_db=self.kmer_db_path,
            fasta=self.fasta,
            crux_comet_params=self.crux_comet_params,
            fasta_fm_index=from_pickle(path=self.fasta_fm_index),
            out_dir=self.hybrid_run_scan_results_dir,
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            min_hybrid_side_len=self.min_hybrid_side_len,
            min_cluster_support=self.min_cluster_support,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
            remove_carbamidomethylated_hybrids=self.remove_carbamidomethylation_hybrids,
            hybrid_decoy_competition=self.hybrid_decoy_competition,
        )


def hybrid_run_on_spectrum(
    spectrum: Spectrum,
    params: Union[HybridRunParams, str, Path],
    fasta_dir: Path,
    on_singularity: bool = False,
    crux_path: Optional[Path] = None,
    num_threads: int = 1,
    delete_hybrids_fasta: bool = True,
):
    # Get params and constants
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    hybrids_fasta = fasta_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    if isinstance(params, (str, Path)):
        params = from_pickle(path=params)
    comet_run = CometRun(
        fasta=hybrids_fasta,
        mzml=spectrum.mzml,
        crux_comet_params=params.crux_comet_params,
        out_dir=params.out_dir,
        decoy_search=2 if params.hybrid_decoy_competition else 0,
        scan_min=spectrum.scan,
        scan_max=spectrum.scan,
        num_threads=num_threads,
    )

    # Form hybrids
    seq_to_hybrids = form_spectrum_hybrids_via_clustering(
        spectrum=spectrum,
        kmer_db=KmerDatabase(db_path=params.kmer_db),
        fasta=Fasta(path=params.fasta),
        fasta_fm_index=params.fasta_fm_index,
        precursor_mz_ppm_tol=params.precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=params.peak_to_ion_ppm_tol,
        min_side_len=params.min_hybrid_side_len,
        min_cluster_support=params.min_cluster_support,
        max_allowed_ion_charge=params.max_allowed_ion_charge,
        remove_carbamidomethylated_hybrids=params.remove_carbamidomethylated_hybrids,
    )
    if len(seq_to_hybrids) == 0:
        logger.info(
            f"No hybrids found for spectrum {spectrum.uid}. Skipping but creating empty Comet outputs..."
        )
        outputs = comet_run.standardized_comet_outputs
        outputs.target.touch()
        if outputs.decoy:
            outputs.decoy.touch()
        return None, comet_run

    # Create hybrids FASTA
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    hybrids_fasta = fasta_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    if params.hybrid_decoy_competition:
        create_hybrids_fasta(
            seq_to_hybrids=seq_to_hybrids,
            output_fasta_path=hybrids_fasta,
            fasta_to_include=params.fasta,
        )
    else:
        create_hybrids_fasta(
            seq_to_hybrids=seq_to_hybrids,
            output_fasta_path=hybrids_fasta,
        )
    process = comet_run.run_comet_and_keep_only_results(
        crux_path=crux_path, on_singularity=on_singularity
    )
    if delete_hybrids_fasta:
        logger.info("Deleting hybrids FASTA")
        os.remove(hybrids_fasta)

    return process, comet_run


def create_hybrids_fasta(
    seq_to_hybrids: Dict[str, List[HybridPeptide]],
    output_fasta_path: Path,
    fasta_to_include: Optional[Path] = None,
    protein_names: Optional[Union[List[str], Path]] = None,
) -> List[Peptide]:
    """
    Writes a list of hybrids (hybrids) to a FASTA file (new_fasta_path).
    If old_fasta is provided, it will also include the proteins from the old FASTA.
    If protein_names is provided, it will only include those proteins from the old FASTA.
    """
    prots = []
    if fasta_to_include is not None:
        if protein_names is not None:
            # Get specific proteins from FASTA by name
            prots = get_proteins_by_name(
                protein_names=protein_names, fasta_path=fasta_to_include
            )
        else:
            # Get all the proteins in the FASTA
            prots = Peptide.from_fasta(fasta_path=fasta_to_include)

    prots.extend(
        HybridPeptide.seq_to_hybrids_map_to_peptides(seq_to_hybrids=seq_to_hybrids)
    )
    Fasta.write_fasta(peptides=prots, path=output_fasta_path)
    return prots


def run_hypedsearch_in_parallel(
    config: Path,
    n_cores: int,
    on_singularity: bool = False,
    crux_path: Optional[Path] = None,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)

    # Get spectra to run Hypedsearch on
    missing_spectra_paths = list(hs_config.missing_hybrid_run_scan_target_txts)
    logger.info(
        f"There are {len(missing_spectra_paths)} spectra missing hybrid Comet outputs."
    )
    if len(missing_spectra_paths) == 0:
        return

    logger.info("Gathering missing spectra")
    missing_spectra = []
    for path in missing_spectra_paths:
        match = re.match(r"^(.+)\.comet\.(\d+)-(\d+)\.(target|decoy)$", path.stem)
        if not match:
            raise ValueError(f"Filename does not match expected pattern: {path}")
        name, start_str, end_str, _ = match.groups()
        assert start_str == end_str
        missing_spectra.append(
            hs_config.spectrum_uid_to_spectrum[
                Spectrum.get_uid(sample=name, scan=int(start_str))
            ]
        )
    logger.info("Done gathering missing spectra.")

    # Pickle shared params
    logger.info("Pickling shared params...")
    params_path = hs_config.parent_output_dir / f"{hs_config.name}_shared_params.pkl"
    to_pickle(obj=hs_config.hybrid_run_params, path=params_path)
    logger.info("Done pickling shared params.")

    # Run Hypedsearch
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create FASTA containing hybrids and run Comet
        logger.info(f"Running Hypedsearch with FASTA dir: {tmp_dir}")
        process_partial = partial(
            hybrid_run_on_spectrum,
            params=params_path,
            on_singularity=on_singularity,
            crux_path=crux_path,
            fasta_dir=Path(tmp_dir),
        )
        with ProcessPoolExecutor(max_workers=n_cores) as ex:
            futures = [
                ex.submit(process_partial, spectrum) for spectrum in missing_spectra
            ]
            # Process results as they complete, catching failures
            results = []
            failed = []
            num_futures = len(futures)
            for idx, future in enumerate(as_completed(futures)):
                try:
                    result = future.result()  # Blocks until THIS task completes
                    results.append(future.result())
                except Exception as e:
                    failed.append(str(e))  # Log failure, continue
                    logger.info(f"Task failed: {e}")

        if len(failed) > 0:
            raise RuntimeError(
                f"{len(failed)}/{num_futures} tasks failed. Errors: {failed}"
            )


@click.command(
    name="combine-comet-txts",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    Given a directory containing the Comet output for each scan, combine them by mzML.
    For each mzML file, two files will be created in out_dir:\n
    1) <mzML name>.comet.target.txt - containing all the target PSMs, and\n
    2) <mzML name>.comet.decoy.txt - containing all the decoy PSMs.
    """,
)
@click.option(
    "--results_dir",
    "-rd",
    type=PathType(),
    required=False,
    help="Path to the directory containing scan results",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=False,
    help="Path to the output directory where the combined results for each mzML will be saved",
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=False,
    help="Path to the Hypedsearch config JSON",
)
def cli_combine_comet_txts(
    results_dir: Optional[Path], out_dir: Optional[Path], config: Optional[Path]
):
    if config is None:
        HypedsearchOutputs.combine_comet_results_by_mzml_and_psm_type(
            folder=results_dir, out_dir=out_dir
        )
    else:
        hs_config = HypedsearchRunConfig.from_json(path=config)
        hs_config.combine_hybrid_run_scan_outputs()
    logger.info("Finished combining Comet scan results")


@click.command(
    name="check-for-missing-scans",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    Check which scans were and were not processed.
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
def cli_check_for_missing_scans(config: Path):
    logger.info(
        f"Checking for missing hybrid scan target txt files for config: {config.name}"
    )
    hs_config = HypedsearchRunConfig.from_json(path=config)
    missing_txts = hs_config.missing_hybrid_run_scan_target_txts
    logger.info(f"Found {len(missing_txts)} missing hybrid scan target txt files")


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_combine_comet_txts)
    cli.add_command(cli_check_for_missing_scans)
    cli()
