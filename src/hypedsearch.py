import logging
import re
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import pandas as pd
from pydantic import BaseModel, field_validator, model_validator
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
    HS_PREFIX,
    HUMAN_PROTEOME,
    HYBRID,
    NATIVE,
    TARGET,
    TRUE_HYBRIDS_PATH,
)
from src.crux import (
    CometConfig,
    CometOutputs,
    Crux,
    get_expected_comet_outputs_for_mzml_to_scans,
)
from src.hybrids_via_clusters import HybridPeptide, form_spectrum_hybrids_via_clustering
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum, create_sample_scan_to_spectrum_map
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.plot_utils import fig_setup, save_fig
from src.protein_abundance import ProteinAbundance
from src.psm import CometPSM
from src.utils import (
    PathType,
    flatten_list_of_lists,
    load_json,
    mass_difference_in_ppm,
    save_dict,
    save_pydantic_objects_to_json,
    setup_logger,
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


def get_matching_output_txt(out_dir: Path, mzml_name: str, psm_type: str) -> Path:
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
                    hybrid_seqs=seq_to_hybrids.keys(),
                    output_fasta_path=hybrids_fasta,
                    fasta_to_include=self.fasta,
                )
            else:
                create_hybrids_fasta(
                    hybrid_seqs=seq_to_hybrids.keys(),
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
            min_cluster_len=self.min_cluster_len,
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
    kmer_db: Optional[Path] = None
    hybrid_decoy_competition: bool = False
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE
    fasta: Path = HUMAN_PROTEOME
    num_peaks: int = 0
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE

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
        # self.plots_dir.mkdir(parents=True, exist_ok=True)
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

    # @property
    # def kmer_db(self) -> Path:
    #     return self.parent_output_dir / f"{self.name}.kmers.db"

    @cached_property
    def spectrum_selector(self) -> SpectrumSelector:
        return SpectrumSelector(max_precursor_charge=self.max_precursor_charge)

    @cached_property
    def hybrid_former(self) -> HybridFormer:
        return HybridFormer(
            kmer_db=self.kmer_db,
            fasta=self.fasta,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            min_cluster_len=self.min_cluster_len,
            min_cluster_support=self.min_cluster_support,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
        )

    @cached_property
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
        return self.native_run_dir / f"{self.name}-{ASSIGN_CONFIDENCE}.txt"

    @property
    def plots_dir(self) -> Path:
        return self.parent_output_dir / "plots"

    @cached_property
    def _mzml_to_scans(self) -> Dict[Path, List[int]]:
        mzml_to_scans = self.mzml_to_scans.copy()
        for mzml, scans in self.mzml_to_scans.items():
            if isinstance(scans, str):
                logger.info
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
            logger.info(f"Native Comet run on MZML {mzml.name} ({idx+1}/{len(mzmls)})")
            outputs.append(
                crux.run_comet(
                    mzml=mzml,
                    fasta=self.hybrid_former.fasta,
                    crux_comet_params=self.psm_scorer.comet_params,
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

    def get_expected_native_run_output_txts(self):
        expected_outputs = self.native_comet_run(dry_run=True)
        results = {
            TARGET: [],
            DECOY: [],
        }
        for output in expected_outputs:
            results[TARGET].append(output.target)
            if output.decoy:
                results[DECOY].append(output.decoy)
        return results

    # def get_expected_hybrid_run_output_txts(self):

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
            prots = prot_cnts.most_common(n=top_n_prots)
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
        )
        return db

    def hybrid_run_on_spectrum(
        self,
        spectrum: Spectrum,
    ) -> Tuple[CometOutputs, CometOutputs, Dict[str, List[HybridPeptide]]]:
        return hybrid_run_on_spectrum(
            spectrum=spectrum,
            out_dir=self.hybrid_run_scan_results_dir,
            hybrid_former=self.hybrid_former,
            psm_scorer=self.psm_scorer,
        )

    def get_output_txts(
        self,
        run_type: Literal[NATIVE, HYBRID],
        psm_type: Literal[TARGET, DECOY, ASSIGN_CONFIDENCE],
    ) -> List:
        if run_type == NATIVE:
            out_dir = self.native_run_dir
        elif run_type == HYBRID:
            out_dir = self.hybrid_run_dir
        else:
            raise ValueError(
                f"Invalid run_type: {run_type}. Allowed: {NATIVE}, {HYBRID}"
            )
        output_txts = []
        for mzml in self.mzml_to_scans.keys():
            mzml_name = Mzml.get_mzml_name(mzml=mzml)
            try:
                output_txts.append(
                    get_matching_output_txt(
                        out_dir=out_dir,
                        mzml_name=mzml_name,
                        psm_type=psm_type,
                    )
                )
            except:
                logger.info(
                    f"Wasn't able to find output txt for MZML={mzml_name}, run_type={run_type}, psm_type={psm_type}"
                )

        return output_txts


# class HypedsearchRunConfig(BaseModel):
#     mzml_to_scans: Dict[Path, Union[Set[int], Literal["all"]]]
#     hybrid_former: HybridFormer
#     psm_scorer: HybridPSMScorer
#     name: str = Field(
#         default_factory=lambda: f"autoGeneratedName_time={datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
#     )
#     native_run_dir: Optional[Path] = None
#     hybrid_run_dir: Optional[Path] = None
#     parent_out_dir: Optional[Path] = None
#     spectrum_selector: SpectrumSelector = SpectrumSelector()
#     spectrum_preprocessor: SpectrumPreprocessor = SpectrumPreprocessor()

#     @field_validator("mzml_to_scans", mode="after")
#     def ensure_mzmls_exist(
#         cls, mzml_to_scans: Dict[Path, Union[List[int], Literal["all"]]]
#     ):
#         for mzml in mzml_to_scans.keys():
#             if not mzml.exists():
#                 raise ValueError(f"MZML file does not exist: {mzml}")
#         return mzml_to_scans

#     @model_validator(mode="after")
#     def post_init(self) -> Self:
#         # Handle native run and hybrid run directories
#         if self.parent_out_dir is None:
#             self.native_run_dir.mkdir(parents=True, exist_ok=True)
#             self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
#         else:
#             self.native_run_dir = self.parent_out_dir / "native_run"
#             self.native_run_dir.mkdir(parents=True, exist_ok=True)
#             self.hybrid_run_dir = self.parent_out_dir / "hybrid_run"
#             self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
#         self.hybrid_run_scan_results_dir.mkdir(parents=True, exist_ok=True)
#         self.native_run_scan_results_dir.mkdir(parents=True, exist_ok=True)
#         return self

#     @cached_property
#     def _mzml_to_scans(self) -> Dict[Path, List[int]]:
#         mzml_to_scans = self.mzml_to_scans.copy()
#         for mzml, scans in self.mzml_to_scans.items():
#             if isinstance(scans, str):
#                 logger.info
#                 mzml_to_scans[mzml] = (
#                     self.spectrum_selector.get_scan_numbers_of_selected_spectra_from_mzml(
#                         mzml=mzml
#                     )
#                 )
#         return mzml_to_scans

#     @cached_property
#     def spectrum_uid_to_spectrum(self) -> Dict[str, Spectrum]:
#         spectrum_uid_to_spectrum = {}
#         for mzml in self.mzml_to_scans.keys():
#             mzml = Mzml(path=mzml)
#             spectrum_uid_to_spectrum.update(mzml.id_to_spectrum)
#         return spectrum_uid_to_spectrum

#     @cached_property
#     def hybrid_run_scan_results_dir(self) -> Path:
#         return self.hybrid_run_dir / "scan_results"

#     @cached_property
#     def native_run_scan_results_dir(self) -> Path:
#         return self.native_run_dir / "scan_results"

#     @classmethod
#     def from_json(cls, path: Union[Path, str]):
#         data = load_json(path=path)
#         return cls(**data)

#     @classmethod
#     def from_data_and_results_dir(
#         cls,
#         data_dir: Union[Path, str],
#         results_dir: Union[Path, str],
#         fasta: Union[Path, str] = HUMAN_PROTEOME,
#         hybrid_native_competition: bool = False,
#     ):
#         data_dir = Path(data_dir)
#         results_dir = Path(results_dir)
#         inputs_dir = results_dir / "inputs"
#         assert inputs_dir.exists(), f"Inputs directory does not exist: {inputs_dir}"
#         comet_params = list(inputs_dir.glob("*comet.params"))
#         assert (
#             len(comet_params) == 1
#         ), f"Expected exactly one comet.params file in {inputs_dir}, found {len(comet_params)}"
#         comet_params = comet_params[0]
#         psm_scorer = {
#             "comet_params": str(comet_params),
#         }
#         if hybrid_native_competition:
#             psm_scorer["fasta"] = str(fasta)
#         hybrid_former = {
#             "kmer_db": str(inputs_dir / "kmers.db"),
#             "fasta": str(fasta),
#         }
#         hs_config = {
#             "name": data_dir.name,
#             "parent_out_dir": str(results_dir),
#             "mzml_to_scans": {
#                 str(mzml): "all" for mzml in list(Path(data_dir).glob("*.mzML"))
#             },
#             "psm_scorer": psm_scorer,
#             "hybrid_former": hybrid_former,
#         }
#         to_json(data=hs_config, path=Path(inputs_dir) / "hs.config.json")

#     def to_json(self, path: Union[str, Path]) -> None:
#         data_as_json_str = self.model_dump_json()
#         to_json(data=json.loads(data_as_json_str), path=path)

#     def native_run(
#         self,
#         method: Literal["direct", "snakemake"] = "direct",
#         native_config: Optional[Path] = None,
#     ):
#         if method == "direct":
#             comet_outputs = self.native_run_via_comet()
#             return
#         elif method == "snakemake":
#             if native_config is None:
#                 native_config = self.native_run_dir / DEFAULT_NATIVE_RUN_CONFIG
#             self.create_native_run_comet_config(
#                 out_path=native_config,
#             )
#             logger.info(
#                 "To run Comet via snakemake, use the following command:\n"
#                 + f"snakemake -s {RUN_COMET_SMK} --configfile {native_config} --scheduler greedy"
#             )

#     def native_run_via_comet(self) -> List[CometOutputs]:
#         crux = Crux()
#         outputs = []
#         mzmls = list(self.mzml_to_scans.keys())
#         for idx, mzml in enumerate(mzmls):
#             logger.info(f"Processing MZML: {mzml.name} ({idx+1}/{len(mzmls)})")
#             outputs.append(
#                 crux.run_comet(
#                     mzml=mzml,
#                     fasta=self.hybrid_former.fasta,
#                     crux_comet_params=self.psm_scorer.comet_params,
#                     decoy_search=2,
#                     out_dir=self.native_run_dir,
#                     file_root=Mzml.get_mzml_name(mzml=mzml),
#                 )
#             )
#         return outputs

#     @staticmethod
#     def assign_confidence_name(prefix: Literal[NATIVE, HYBRID]) -> str:
#         return f"{prefix}-assign-confidence.txt"

#     @property
#     def native_assign_confidence_path(self) -> Path:
#         return self.native_run_dir / self.assign_confidence_name(prefix=NATIVE)

#     @property
#     def hybrid_assign_confidence_path(self) -> Path:
#         return self.hybrid_run_dir / self.assign_confidence_name(prefix=HYBRID)

#     def run_assign_confidence(self, run_type: Literal[NATIVE, HYBRID] = NATIVE):
#         if run_type == NATIVE:
#             out_dir = self.native_run_dir
#             out_path = self.native_assign_confidence_path
#         elif run_type == HYBRID:
#             out_dir = self.hybrid_run_dir
#             out_path = self.hybrid_assign_confidence_path
#         else:
#             raise ValueError(
#                 f"Invalid run_type: {run_type}. Allowed: {NATIVE}, {HYBRID}"
#             )
#         target_txts = list(out_dir.glob(f"*.{TARGET}.txt"))
#         decoy_txts = list(out_dir.glob(f"*.{DECOY}.txt"))
#         crux = Crux()
#         if len(decoy_txts) == 0:
#             logger.info(f"No decoy txt files found in {out_dir}. Skipping.")
#             return
#         crux.run_assign_confidence(
#             target_txts=target_txts,
#             out_path=out_path,
#         )

#     @cached_property
#     def spectra(self) -> List[Spectrum]:
#         spectra = []
#         for mzml, scans in self._mzml_to_scans.items():
#             mzml_spectra = Spectrum.parse_ms2_from_mzml(mzml=mzml)
#             spectra_from_mzml = [
#                 spectrum for spectrum in mzml_spectra if spectrum.scan in scans
#             ]
#             spectra.extend(spectra_from_mzml)
#         return spectra

#     def native_run_on_spectrum(
#         self, spectrum: Spectrum, num_threads_4_comet: int = DEFAULT_NUM_COMET_THREADS
#     ) -> CometOutputs:
#         comet_config = self.create_native_run_comet_config(
#             num_threads_4_comet=num_threads_4_comet, decoy_search=2
#         )
#         outputs = native_run_on_spectrum(spectrum=spectrum, comet_config=comet_config)
#         return outputs

#     def run_hypedsearch_on_spectrum(
#         self, spectrum: Spectrum, num_threads_4_comet: int = DEFAULT_NUM_COMET_THREADS
#     ) -> Tuple[CometOutputs, CometOutputs, Dict[str, List[HybridPeptide]]]:
#         native_outputs = self.native_run_on_spectrum(
#             spectrum=spectrum, num_threads_4_comet=num_threads_4_comet
#         )
#         seq_to_hybrids, hybrid_outputs = hybrid_run_on_spectrum(
#             spectrum=spectrum,
#             out_dir=self.hybrid_run_scan_results_dir,
#             spectrum_preprocessor=self.spectrum_preprocessor,
#             hybrid_former=self.hybrid_former,
#             psm_scorer=self.psm_scorer,
#         )

#         return native_outputs, hybrid_outputs, seq_to_hybrids

#     @property
#     def expected_native_scan_target_outputs(self) -> List[str]:
#         return get_expected_comet_outputs_for_mzml_to_scans(
#             mzml_to_scans=self._mzml_to_scans,
#             out_dir=self.native_run_scan_results_dir,
#             decoy_search=2,  # this value doesn't matter
#             psm_type=TARGET,
#         )

#     @property
#     def missing_native_target_outputs(self) -> List[str]:
#         actual = set([txt for txt in self.native_run_scan_results_dir.glob("*.txt")])
#         expected = set(self.expected_native_scan_target_outputs)
#         return [str(txt) for txt in (expected - actual)]

#     @property
#     def expected_hybrid_scan_target_outputs(self) -> List[str]:
#         return self.psm_scorer.expected_hybrid_target_outputs(
#             mzml_to_scans=self._mzml_to_scans,
#             out_dir=self.hybrid_run_scan_results_dir,
#         )

#     @property
#     def missing_hybrid_target_outputs(self) -> List[str]:
#         actual = set([txt for txt in self.hybrid_run_scan_results_dir.glob("*.txt")])
#         expected = set(self.expected_hybrid_scan_target_outputs)
#         return [str(txt) for txt in (expected - actual)]

#     def print_missing_scan_info(self) -> Tuple[List[str], List[str]]:
#         missing_hybrid_targets = self.missing_hybrid_target_outputs
#         missing_native_targets = self.missing_native_target_outputs
#         print(
#             f"There are\n- {len(missing_hybrid_targets)} missing hybrid scan (target) outputs and\n- {len(missing_native_targets)} missing native scan (target) outputs."
#         )
#         return missing_native_targets, missing_hybrid_targets

#     def combine_comet_scan_results(self):
#         num_missing_natives = len(self.missing_native_target_outputs)
#         num_missin_hybrids = len(self.missing_hybrid_target_outputs)
#         assert (num_missing_natives == 0) and (num_missin_hybrids == 0), (
#             f"There are {num_missing_natives} missing native scan results and {num_missin_hybrids} missing hybrids! "
#             "Aborting because you probably want to have every spectrum's results before combining them by MZML."
#         )
#         # Combine native results
#         logger.info("Combining native scan results...")
#         combine_comet_scan_results(
#             scan_results_dir=self.native_run_scan_results_dir,
#             out_dir=self.native_run_dir,
#         )

#         # Combine hybrid results
#         logger.info("Combining hybrid scan results...")
#         combine_comet_scan_results(
#             scan_results_dir=self.hybrid_run_scan_results_dir,
#             out_dir=self.hybrid_run_dir,
#         )

# def get_output_txts(
#     self,
#     run_type: Literal[NATIVE, HYBRID],
#     psm_type: Literal[TARGET, DECOY, ASSIGN_CONFIDENCE],
# ) -> List:
#     def file_to_load(run_type, psm_type):
#         if psm_type == TARGET:
#             return f"*.comet.{TARGET}.txt"
#         elif psm_type == DECOY:
#             return f"*.comet.{DECOY}.txt"
#         elif psm_type == ASSIGN_CONFIDENCE:
#             return self.assign_confidence_name(prefix=run_type)
#         else:
#             raise ValueError(
#                 f"Invalid psm_type: {psm_type}. Allowed: {TARGET}, {DECOY}, {ASSIGN_CONFIDENCE}"
#             )

#     if run_type == NATIVE:
#         out_dir = self.native_run_dir
#     elif run_type == HYBRID:
#         out_dir = self.hybrid_run_dir
#     else:
#         raise ValueError(
#             f"Invalid run_type: {run_type}. Allowed: {NATIVE}, {HYBRID}"
#         )
#     output_txts = []
#     for mzml in self.mzml_to_scans.keys():
#         mzml_name = Mzml.get_mzml_name(mzml=mzml)
#         try:
#             output_txts.append(
#                 get_matching_output_txt(
#                     out_dir=out_dir,
#                     mzml_name=mzml_name,
#                     psm_type=psm_type,
#                 )
#             )
#         except:
#             logger.info(
#                 f"Wasn't able to find output txt for MZML={mzml_name}, run_type={run_type}, psm_type={psm_type}"
#             )

#     return output_txts

#     @property
#     def _results_dir(self) -> Path:
#         results_dir = self.parent_out_dir / f"{self.name}/{DEFAULT_RESULTS_DIR_NAME}"
#         results_dir.mkdir(parents=True, exist_ok=True)
#         return results_dir


def native_run_on_spectrum(
    spectrum: Spectrum, comet_config: CometConfig
) -> CometOutputs:
    return comet_config.run_comet_on_mzml(
        mzml=spectrum.mzml, scan_min=spectrum.scan, scan_max=spectrum.scan
    )


def hybrid_run_on_spectrum(
    spectrum: Spectrum,
    out_dir: Path,
    hybrid_former: HybridFormer,
    psm_scorer: HybridPSMScorer,
) -> Tuple[Dict[str, List[HybridPeptide]], List[CometOutputs]]:
    seq_to_hybrids = hybrid_former.form_hybrids(spectrum=spectrum)
    comet_outputs = psm_scorer.score_hybrids(
        seq_to_hybrids=seq_to_hybrids,
        spectrum=spectrum,
        out_dir=out_dir,
    )
    return seq_to_hybrids, comet_outputs


def hybrid_fasta_name(hybrid_seq: str) -> str:
    return f"{HS_PREFIX}{hybrid_seq}"


def create_hybrids_fasta(
    hybrid_seqs: List[str],
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

    for _, hybrid_seq in enumerate(hybrid_seqs):
        new_peptide = Peptide(
            seq=hybrid_seq,
            name=hybrid_fasta_name(hybrid_seq=hybrid_seq),
        )

        prots.append(new_peptide)

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


# @click.command(
#     name="create-native-run-snakemake-config",
#     context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
#     help="""
#     Create the native run snakemake config
#     """,
# )
# @click.option(
#     "--config",
#     "-c",
#     type=PathType(),
#     required=True,
#     help="Path to the Hypedsearch config JSON",
# )
# @click.option(
#     "--out_path",
#     "-o",
#     type=PathType(),
#     required=False,
#     help=f"Path to where the native run snakemake config will be saved. Default is <native_run_dir>/{DEFAULT_NATIVE_RUN_CONFIG}",
# )
# @click.option(
#     "--num_threads",
#     "-n",
#     type=int,
#     required=False,
#     default=DEFAULT_NUM_COMET_THREADS,
#     show_default=True,
#     help="Number of threads to use for Comet in the native run",
# )
# def cli_create_native_run_snakemake_config(
#     config: Path, out_path: Optional[Path], num_threads: int
# ):
#     hs_config = HypedsearchRunConfig.from_json(path=config)
#     if out_path is None:
#         out_path = hs_config.native_run_dir / DEFAULT_NATIVE_RUN_CONFIG
#     hs_config.create_native_run_comet_config(
#         out_path=out_path, num_threads_4_comet=num_threads
#     )


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


@click.command(
    name="run-hypedsearch",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    Run HypedSearch on a single scan from an mzML file.
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to Hypedsearch JSON config",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=False,
    help="Path to MZML",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    required=False,
    help="Scan in MZML on which to run HypedSearch",
)
def cli_run_hypedsearch(mzml: Optional[Path], scan: Optional[int], config: Path):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    # Run on one spectrum
    if scan is not None:
        mzml_obj = Mzml(path=mzml)
        return hs_config.hybrid_run_on_spectrum(
            spectrum=mzml_obj.get_spectrum(scan=scan)
        )
    # Run on all spectra specified in config
    for mzml, scans in hs_config._mzml_to_scans.items():
        mzml_obj = Mzml(path=mzml)
        for scan in scans:
            _ = hs_config.hybrid_run_on_spectrum(
                spectrum=mzml_obj.get_spectrum(scan=scan)
            )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_combine_comet_txts)
    cli.add_command(cli_check_for_missing_scans)
    cli.add_command(cli_run_hypedsearch)
    # cli.add_command(cli_create_native_run_snakemake_config)
    cli()
