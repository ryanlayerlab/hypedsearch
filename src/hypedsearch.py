import logging
import multiprocessing as mp
import re
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from functools import cached_property, partial
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import pandas as pd
from fm_index import MultiFMIndex
from pydantic import BaseModel, ConfigDict, field_validator, model_validator
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
    SHARED_PARAMS,
    TARGET,
    TRUE_HYBRIDS_PATH,
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
from src.plot_utils import fig_setup, save_fig
from src.psm import CometPSM, ProteinAbundance
from src.utils import (
    PathType,
    flatten_list_of_lists,
    from_pickle,
    load_json,
    log_params,
    mass_difference_in_ppm,
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
    def get_combined_comet_output_name(mzml_name: str, psm_type: str) -> str:
        return f"{mzml_name}.comet.0-0.{psm_type}.txt"

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


class HybridFormer(BaseModel):
    kmer_db: Path
    fasta: Path
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE

    def form_hybrids(self, spectrum: Spectrum) -> Dict[str, List[HybridPeptide]]:
        seq_to_hybrids = form_spectrum_hybrids_via_clustering(
            spectrum=spectrum,
            kmer_db=KmerDatabase(db_path=self.kmer_db),
            fasta=Fasta(path=self.fasta),
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            min_side_len=self.min_cluster_len,
            min_cluster_support=self.min_cluster_support,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
        )
        return seq_to_hybrids


class SpectrumSelector(BaseModel):
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE

    def select_spectra(self, spectra: List[Spectrum]) -> List[Spectrum]:
        return [
            spectrum
            for spectrum in spectra
            if spectrum.precursor_charge <= self.max_precursor_charge
        ]

    def get_scan_numbers_of_selected_spectra_from_mzml(self, mzml: Path) -> Set[int]:
        logger.info(f"Getting spectra from {mzml.name}")
        spectra = self.select_spectra(spectra=Spectrum.parse_ms2_from_mzml(mzml=mzml))
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


class HypedsearchRunConfig(BaseModel):
    mzml_to_scans: Dict[Path, Union[Set[int], Literal["all"]]]
    parent_output_dir: Path
    crux_comet_params: Path
    name: str
    fasta_fm_index: Path
    kmer_db: Optional[Path] = None
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

    @field_validator("mzml_to_scans", mode="after")
    def ensure_mzmls_exist(
        cls, mzml_to_scans: Dict[Path, Union[List[int], Literal["all"]]]
    ):
        for mzml in mzml_to_scans.keys():
            if not mzml.exists():
                raise ValueError(f"MZML file does not exist: {mzml}")
        return mzml_to_scans

    @field_validator("crux_comet_params", mode="after")
    def ensure_paths_exist(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"Path does not exist: {v}")
        return v

    @model_validator(mode="after")
    def post_init(self) -> Self:
        # Set kmer database path
        if self.kmer_db is None:
            self.kmer_db = self.parent_output_dir / f"{self.name}.kmers.db"

        # Ensure directories exist
        self.native_run_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_run_scan_results_dir.mkdir(parents=True, exist_ok=True)

        return self

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
    def spectrum_uid_to_spectrum(self) -> Dict[str, Spectrum]:
        spectrum_uid_to_spectrum = {}
        for mzml in self.mzml_to_scans.keys():
            mzml = Mzml(path=mzml)
            spectrum_uid_to_spectrum.update(mzml.id_to_spectrum)
        return spectrum_uid_to_spectrum

    @property
    def spectra(self) -> List[Spectrum]:
        return list(self.spectrum_uid_to_spectrum.values())

    @classmethod
    def from_json(cls, path: Union[Path, str]):
        data = load_json(path=path)
        Path(data["parent_output_dir"]).mkdir(exist_ok=True, parents=True)
        return cls(**data)

    @cached_property
    def native_assign_confidence_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.native_assign_confidence_path)

    @property
    def native_assign_confidence_path(self) -> Path:
        return self.native_run_dir / HypedsearchOutputs.get_expected_output_txt_name(
            name=self.name, psm_type=ASSIGN_CONFIDENCE
        )

    @property
    def plots_dir(self) -> Path:
        return self.parent_output_dir / "plots"

    @cached_property
    def _mzml_to_scans(self) -> Dict[Path, List[int]]:
        mzml_to_scans = self.mzml_to_scans.copy()
        for mzml, scans in self.mzml_to_scans.items():
            if isinstance(scans, str):
                mzml_to_scans[mzml] = (
                    self.spectrum_selector.get_scan_numbers_of_selected_spectra_from_mzml(
                        mzml=mzml
                    )
                )
        return mzml_to_scans

    @property
    def expected_hybrid_run_scan_target_txts(self) -> Set[Path]:
        txts = self.psm_scorer.expected_hybrid_target_outputs(
            mzml_to_scans=self._mzml_to_scans,
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
        return self.model_dump(mode="json")

    def save(
        self,
        path: Union[Path, str],
    ):
        """Save config as a JSON"""
        save_dict(data=self.to_dict(), path=path)

    def native_comet_run(self, dry_run: bool = False) -> List[CometOutputs]:
        crux = Crux()
        outputs = []
        mzmls = list(self.mzml_to_scans.keys())
        for idx, mzml in enumerate(mzmls):
            if not dry_run:
                logger.info(
                    f"Native Comet run on MZML {mzml.name} ({idx+1}/{len(mzmls)})"
                )
            outputs.append(
                crux.run_comet(
                    mzml=mzml,
                    fasta=self.fasta,
                    crux_comet_params=self.crux_comet_params,
                    decoy_search=2,
                    out_dir=self.native_run_dir,
                    file_root=Mzml.get_mzml_name(mzml=mzml),
                    dry_run=dry_run,
                )
            )
        return outputs

    def native_assign_confidence(self):
        native_outputs = self.get_expected_native_run_output_txts()
        Crux().run_assign_confidence(
            target_txts=native_outputs[TARGET],
            out_path=self.native_assign_confidence_path,
        )

    def get_expected_native_run_output_txts(self) -> Dict[str, List[Path]]:
        expected_outputs = self.native_comet_run(dry_run=True)
        results = {
            TARGET: [],
            DECOY: [],
            ASSIGN_CONFIDENCE: self.native_run_dir
            / HypedsearchOutputs.get_expected_output_txt_name(
                name=self.name, psm_type=ASSIGN_CONFIDENCE
            ),
        }
        for output in expected_outputs:
            results[TARGET].append(output.target)
            if output.decoy:
                results[DECOY].append(output.decoy)
        return results

    def get_expected_hybrid_run_output_txts(self) -> List[Path]:
        return [
            self.hybrid_run_dir
            / HypedsearchOutputs.get_expected_output_txt_name(
                name=mzml_name, psm_type=TARGET, scan=0
            )
            for mzml_name in self.mzml_names
        ]

    def get_protein_abundance(
        self,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
        top_n_prots: int = 20,
        plot: bool = False,
    ) -> ProteinAbundance:
        prot_ab = ProteinAbundance.from_comet_psms(
            psms=self.native_assign_confidence_psms, q_threshold=q_threshold
        )
        if plot:
            fig, axs = fig_setup(h=8)
            ax = axs[0]
            prot_ab.plot(top_n_prots=top_n_prots, ax=ax)
            save_fig(
                fig=fig,
                title=f"{self.name}\nTop {top_n_prots} proteins by number of PSMs w/ q<={q_threshold}",
                path=self.plots_dir / f"{self.name}.protAb.png",
            )
        return prot_ab

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
            db_path=self.kmer_db,
            proteins=fasta.get_proteins_by_name(names=prots),
            min_k=min_k,
            max_k=max_k,
            overwrite=True,
        )
        return db

    def hybrid_run_on_spectrum(
        self,
        spectrum: Spectrum,
    ) -> Tuple[Dict[str, List[HybridPeptide]], List[CometOutputs]]:
        seq_to_hybrids = form_spectrum_hybrids_via_clustering(
            spectrum=spectrum,
            kmer_db=KmerDatabase(db_path=self.kmer_db),
            fasta=Fasta(path=self.fasta),
            fasta_fm_index=Fasta2MFMIndex.load(path=self.fasta_fm_index),
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            min_side_len=self.min_hybrid_side_len,
            min_cluster_support=self.min_cluster_support,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
            remove_carbamidomethylated_hybrids=self.remove_carbamidomethylation_hybrids,
        )
        comet_outputs = self.psm_scorer.score_hybrids(
            seq_to_hybrids=seq_to_hybrids,
            spectrum=spectrum,
            out_dir=self.hybrid_run_scan_results_dir,
        )
        return seq_to_hybrids, comet_outputs


@dataclass
class HybridRunParams:
    kmer_db: Path
    fasta: Path
    crux_comet_params: Path
    fasta_fm_index: MultiFMIndex
    out_dir: Path
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    min_hybrid_side_len: int = DEFAULT_MIN_SIDE_LEN
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE
    remove_carbamidomethylated_hybrids: bool = True
    hybrid_decoy_competition: bool = False


def hybrid_run_on_spectrum(
    spectrum: Spectrum,
    params: Union[HybridRunParams, str, Path],
    on_singularity: bool = False,
    crux_path: Optional[Path] = None,
    num_threads: int = 1,
):
    if isinstance(params, (str, Path)):
        params = from_pickle(path=params)
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
    # Create hybrids FASTA
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create FASTA containing hybrids and run Comet
        tmp_path = Path(tmp_dir)
        mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
        hybrids_fasta = tmp_path / f"{mzml_name}.{spectrum.scan}.fasta"
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
        # Run Comet
        if params.hybrid_decoy_competition:
            decoy_search = 2
        else:
            decoy_search = 0
        run = CometRun(
            fasta=hybrids_fasta,
            mzml=spectrum.mzml,
            crux_comet_params=params.crux_comet_params,
            out_dir=params.out_dir,
            decoy_search=decoy_search,
            scan_min=spectrum.scan,
            scan_max=spectrum.scan,
            num_threads=num_threads,
        )
        process = run.run_comet_and_keep_only_results(
            crux_path=crux_path, on_singularity=on_singularity
        )
    return process, run


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
        HypedsearchOutputs.combine_comet_results_by_mzml_and_psm_type(
            folder=hs_config.hybrid_run_scan_results_dir
        )
    logger.info("Finished combining Comet scan results")


def run_hypedsearch_in_parallel(
    config: Path,
    n_cores: int,
    on_singularity: bool = False,
    crux_path: Optional[Path] = None,
):
    # # Disable threading in child processes (critical for C/C++ libs)
    # mp.set_start_method(
    #     "spawn", force=True
    # )  # Fresh Python interp, no inherited state [web:13][cite:28]

    hs_config = HypedsearchRunConfig.from_json(path=config)
    # Pickle shared params
    logger.info("Pickling shared params...")
    params = HybridRunParams(
        kmer_db=hs_config.kmer_db,
        fasta=hs_config.fasta,
        fasta_fm_index=from_pickle(path=hs_config.fasta_fm_index),
        crux_comet_params=hs_config.crux_comet_params,
        out_dir=hs_config.hybrid_run_scan_results_dir,
    )
    params_path = hs_config.parent_output_dir / f"{hs_config.name}_shared_params.pkl"
    to_pickle(obj=params, path=params_path)

    # Get spectra to run Hypedsearch on
    missing_spectra_paths = list(hs_config.missing_hybrid_run_scan_target_txts)
    logger.info(
        f"There are {len(missing_spectra_paths)} spectra missing hybrid Comet outputs. So I'll run Hypedsearch on those spectra."
    )
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

    # Run Hypedsearch
    logger.info("Running Hypedsearch")
    process_partial = partial(
        hybrid_run_on_spectrum,
        params=params_path,
        on_singularity=on_singularity,
        crux_path=crux_path,
    )
    with ProcessPoolExecutor(max_workers=n_cores) as ex:
        # futures = []
        # for i in range(0, len(spectra), chunksize):
        #     batch = spectra[i : i + chunksize]
        #     futures.extend(ex.submit(process_partial, s) for s in batch)
        # results = list(ex.map(process_partial, spectra, chunksize=chunksize))
        futures = [ex.submit(process_partial, spectrum) for spectrum in missing_spectra]

        # Process results as they complete, catching failures
        results = []
        failed = []
        for future in as_completed(futures):
            try:
                result = future.result()  # Blocks until THIS task completes
                results.append(result)
            except Exception as e:
                failed.append(str(e))  # Log failure, continue
                logger.info(f"Task failed: {e}")

    logger.info(f"Completed: {len(results)}, Failed: {len(failed)}")


@click.command(
    name="run-in-parallel",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    required=True,
    help="",
)
@click.option(
    "--on_singularity",
    "-os",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@log_params
def cli_run_in_parallel(config: Path, n_cores: int, on_singularity: bool):
    if on_singularity:
        crux_path = None
    else:
        crux_path = MAC_CRUX_EXECUTABLE
    run_hypedsearch_in_parallel(
        config=config,
        n_cores=n_cores,
        crux_path=crux_path,
        on_singularity=on_singularity,
    )


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
    cli.add_command(cli_run_in_parallel)
    # cli.add_command(cli_create_native_run_snakemake_config)
    cli()
