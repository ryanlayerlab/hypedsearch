import logging
import os
import re
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from functools import cached_property, partial
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pydantic import BaseModel, model_validator
from typing_extensions import Self

from src.constants import (
    ALL,
    ASSIGN_CONFIDENCE,
    COMET_PROTEIN_SEPARATOR,
    DECOY,
    DEFAULT_FPR,
    DEFAULT_MAX_ALLOWED_ION_CHARGE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MAX_PRECURSOR_CHARGE,
    DEFAULT_MIN_CLUSTER_SUPPORT,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_Q_RANGE,
    DEFAULT_Q_THRESHOLD,
    DEFAULT_SCORE_CHANGE_RANGE,
    MAC_CRUX_EXECUTABLE,
    PLAIN_PEPTIDE,
    PROTEIN,
    TARGET,
    TRUE_HYBRIDS_PATH,
)
from src.crux import CometOutputs, CometRun, Crux
from src.hybrids_via_clusters import HybridPeptide, form_spectrum_hybrids_via_clustering
from src.kmer_database import KmerDatabase
from src.mass_spectra import (
    Mzml,
    Spectrum,
    create_sample_scan_to_spectrum_map,
    organize_by_spectrum_uid,
    plot_spectra_histograms,
)
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.psm import (
    CometPSM,
    CometRunAnalysis,
    ProteinAbundance,
    create_xcorr_dists_plot,
)
from src.utils import (
    CmdLineResult,
    check_if_file_is_empty,
    flatten_list_of_lists,
    from_pickle,
    load_json,
    log_params,
    log_time,
    mass_difference_in_ppm,
    path_aware_dict_factory,
    save_dict,
    save_pydantic_objects_to_json,
    setup_logger,
    to_json,
)

logger = logging.getLogger(__name__)


def find_possible_hybrids_for_seq(
    seq: str, kmer_to_proteins_map: Dict[str, List[str]], min_side_len: int
) -> List[HybridPeptide]:
    """Find the possible hybrid peptides formed from the given k-mers that would share the sequence `seq`.

    Args:
        seq: Sequence to find hybrids for
        kmer_to_proteins_map: Mapping from k-mers to the proteins they appear in.
        min_side_len: Require hybrid peptides to have at least this many amino acids on each side.

    Returns:
        A list of hybrids that formed from `kmer_to_proteins_map` that share the sequence `seq`
            and have at least `min_side_len` amino acids on each side.
    """
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


def get_seq_to_hybrids_map(
    seqs: Union[Set[str], List[str]],
    kmer_db_path: Path,
    min_side_len: int = DEFAULT_MIN_SIDE_LEN,
    remove_carbamidomethylation: bool = True,
) -> Dict[str, List[HybridPeptide]]:
    """Given a set of sequences and a k-mer database, find the hybrids from the kmer-

    Args:
        seqs: A set of sequences to find hybrids for.
        kmer_db_path: Path to k-mer database.
        min_side_len: Require hybrid peptides to have at least this many amino acids on each side.
            Defaults to DEFAULT_MIN_SIDE_LEN.
        remove_carbamidomethylation: Whether or not to remove hybrids that show evidence
            of carbamidomethylation. Defaults to True.

    Returns:
        A mapping from each sequence `seq` in `seqs` to the list of hybrids as
        `HybridPeptide` objects that could form `seq`.
    """
    kmer_to_proteins_map = KmerDatabase(
        db_path=kmer_db_path
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


class HybridRunParams(BaseModel):
    """Class for the parameters that define a HypedSearch run.

    Args:
        kmer_db_path: Path to k-mer database.
        fasta: Path to FASTA file.
        fasta_fm_index: Path to the FM index for the FASTA file.
        crux_comet_params: Path to the parameter file to pass [`crux comet`](https://crux.ms/commands/comet.html)
            for both the native and hybrid Comet runs.
        hybrid_decoy_competition: Whether or not to perform a target-decoy competition in
            the hybrid Comet run.
        max_precursor_chrage: Only spectra with precursor charge less than this will be
            analyzed.
        peak_to_ion_ppm_tol: When forming hybrids and matching theoretical product ions
            to spectrum peaks, require that the PPM difference between the observed peak
            and the theoretical ion is less than this tolerance.
        precursor_mz_ppm_tol: When forming hybrids, require that the PPM difference between
            the spectrum precursor m/z and the theoretical hybrid precursor m/z is less
            than this tolerance.
        min_hybrid_side_len: When forming hybrids, require that each side of the hybrid
            has at least this many amino acids.
        min_cluster_support: When forming hybrids via b- and y-clusters, require that b-
            and y-clusters have at least this many supporting ions to be considered.
        max_allowed_ion_charge: When forming hybrids, only consider b- and y-ions with
            charge less than or equal to this.
        remove_carbamidomethylation_hybrids: Whether or not to remove hybrids that show
            evidence of carbamidomethylation.
    """

    kmer_db_path: Path
    fasta: Path
    fasta_fm_index: Path
    crux_comet_params: Path
    hybrid_decoy_competition: bool = False
    max_precursor_charge: int = DEFAULT_MAX_PRECURSOR_CHARGE
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    min_hybrid_side_len: int = DEFAULT_MIN_SIDE_LEN
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE
    remove_carbamidomethylation_hybrids: bool = True

    @property
    def kmer_db(self) -> KmerDatabase:
        return KmerDatabase(db_path=self.kmer_db_path)

    @classmethod
    def load(cls, path: str | Path) -> "HybridRunParams":
        """Load and return a `HybridRunParams` object from a JSON file."""
        data = load_json(path=path)
        return cls(**data)

    def save(self, path: str | Path):
        """Save the HybridRunParams object as a JSON file."""
        to_json(
            data=self.model_dump(mode="json"),
            path=path,
        )

    def form_hybrids(
        self,
        spectrum: Spectrum,
    ):
        return form_spectrum_hybrids_via_clustering(
            spectrum=spectrum,
            kmer_db=self.kmer_db,
            fasta=Fasta(path=self.fasta),
            fasta_fm_index=from_pickle(self.fasta_fm_index),
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            min_side_len=self.min_hybrid_side_len,
            min_cluster_support=self.min_cluster_support,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
            remove_carbamidomethylated_hybrids=self.remove_carbamidomethylation_hybrids,
        )


@dataclass
class HypedsearchRunConfig:
    """Helper class for managing a HypedSearch run. Contains properties and methods for
    getting spectra, running Comet natively, running HypedSearch, and analyzing results.

    Args:
        name: Name of the HypedSearch run/configuration. Used for naming output files and directories.
        mzml_to_scans: Mapping from MZML file paths to the set of scan numbers to run
            HypedSearch on for that MZML.
        parent_output_dir: Parent directory under which output and results directories
            will be created under `<parent_output_dir>/<name>`
        hybrid_run_params: Parameters for the HypedSearch run.
    """

    name: str
    mzml_to_scans: Dict[Path, Union[Set[int], Literal[ALL]]]
    parent_output_dir: Path
    hybrid_run_params: HybridRunParams

    def __post_init__(self):
        # Set any strings that are supposed to be Path objects to Path objects
        self.mzml_to_scans = {
            Path(mzml): scans for mzml, scans in self.mzml_to_scans.items()
        }
        self.parent_output_dir = Path(self.parent_output_dir)

        return self

    @property
    def kmer_db(self) -> KmerDatabase:
        return self.hybrid_run_params.kmer_db

    @property
    def mzml_names(self) -> List[str]:
        return [Mzml(path=mzml).name for mzml in self.mzml_to_scans.keys()]

    # Spectra properties
    @cached_property
    def _mzml_to_spectra(self) -> Dict[Path, List[Spectrum]]:
        logger.info(f"Getting spectra for config {self.name}")
        mzml_to_spectra = {}
        for mzml, scans in self.mzml_to_scans.items():
            all_spectra = Spectrum.parse_ms2_from_mzml(mzml=mzml)
            if scans == ALL:
                mzml_to_spectra[mzml] = all_spectra
            else:
                # Filter to only spectra with scan numbers in scans
                mzml_to_spectra[mzml] = [
                    spectrum for spectrum in all_spectra if spectrum.scan in scans
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
    def uid_to_spectrum(self) -> Dict[str, Spectrum]:
        return self.spectrum_uid_to_spectrum

    @property
    def spectra(self) -> List[Spectrum]:
        return list(self.spectrum_uid_to_spectrum.values())

    @cached_property
    def spectra_df(self):
        return Spectrum.to_df(
            spectra=self.spectra,
        )

    @cached_property
    def _mzml_to_scan_nums(self) -> Dict[Path, List[int]]:
        mzml_to_scan_nums = {}
        for mzml, spectra in self._mzml_to_spectra.items():
            mzml_to_scan_nums[mzml] = [spectrum.scan for spectrum in spectra]
        return mzml_to_scan_nums

    @cached_property
    def spectrum_uids_meeting_selection_criteria(self) -> Set[str]:
        uids = set()
        for uid, spectrum in self.uid_to_spectrum.items():
            if spectrum.precursor_charge <= self.hybrid_run_params.max_precursor_charge:
                uids.add(uid)
        return uids

    @property
    def name_dir(self) -> Path:
        d = self.parent_output_dir / f"{self.name}"
        d.mkdir(parents=True, exist_ok=True)
        return d

    # Native run properties
    @cached_property
    def native_run_dir(self) -> Path:
        d = self.name_dir / "native_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def expected_native_comet_outputs(self) -> List[Path]:
        expected_outputs = []
        for mzml_path in self.mzml_to_scans.keys():
            comet_run = CometRun(
                fasta=self.hybrid_run_params.fasta,
                mzml=mzml_path,
                crux_comet_params=self.hybrid_run_params.crux_comet_params,
                decoy_search=2,
                out_dir=self.native_run_dir,
                file_root=Mzml.get_mzml_name(mzml=mzml_path),
            )
            expected_outputs.append(comet_run.standardized_comet_outputs)
        return expected_outputs

    @cached_property
    def native_target_psms(self) -> List[CometPSM]:
        logger.info("Getting native target PSMs...")
        txts = []
        for output in self.expected_native_comet_outputs:
            assert (
                output.target.exists()
            ), f"Expected target txt does not exist: {output.target}"
            txts.append(output.target)
        psms = CometPSM.from_txts(txts=txts)
        return psms

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
    def native_results_dir(self) -> Path:
        d = self.parent_output_dir / f"results/{self.name}_native_results"
        d.mkdir(exist_ok=True, parents=True)
        return d

    @property
    def native_assign_confidence_path(self) -> Path:
        return self.native_run_dir / self.get_expected_output_txt_name(
            name=self.name, psm_type=ASSIGN_CONFIDENCE
        )

    @cached_property
    def native_assign_confidence_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.native_assign_confidence_path)

    @cached_property
    def native_comet_run(self):
        return CometRunAnalysis(
            targets=self.native_target_psms,
            decoys=self.native_decoy_psms,
            assign_conf=self.native_assign_confidence_psms,
        )

    # Hybrid run properties
    @property
    def expected_hybrid_txts(self):
        return [
            Path(self.hybrid_run_dir)
            / self.get_combined_comet_output_name(mzml_name=mzml_name, psm_type=TARGET)
            for mzml_name in self.mzml_names
        ]

    @cached_property
    def hybrid_run_dir(self) -> Path:
        d = self.name_dir / "hybrid_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @cached_property
    def hybrid_run_scan_results_dir(self) -> Path:
        d = self.hybrid_run_dir / "scan_results"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @cached_property
    def hybrid_target_psms(self):
        logger.info("Getting hybrid target PSMs...")
        txts = []
        for txt in self.expected_hybrid_txts:
            assert txt.exists(), f"Expected hybrid target txt does not exist: {txt}"
            txts.append(txt)
        psms = CometPSM.from_txts(txts=txts)
        return psms

    @property
    def expected_hybrid_run_spectrum_target_txts(self) -> Set[Path]:

        return self.get_expected_hybrid_comet_outputs_when_running_on_all()

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
            self.expected_hybrid_run_spectrum_target_txts
            - self.existing_hybrid_run_scan_target_txts
        )

    @cached_property
    def hybrid_comet_run(self):
        return CometRunAnalysis(
            targets=self.hybrid_target_psms,
            assign_conf=self.native_assign_confidence_psms,
        )

    # Class methods
    @classmethod
    def from_json(cls, path: Union[Path, str]):
        data = load_json(path=path)
        Path(data["parent_output_dir"]).mkdir(exist_ok=True, parents=True)
        if isinstance(data["hybrid_run_params"], dict):
            data["hybrid_run_params"] = HybridRunParams(**data["hybrid_run_params"])
        else:
            data["hybrid_run_params"] = HybridRunParams.load(
                path=data["hybrid_run_params"]
            )
        return cls(**data)

    # Instance methods
    def get_expected_hybrid_comet_outputs_when_running_on_all(self) -> List[Path]:
        expected_outputs = set()
        for psm in self.native_target_psms:
            if psm.uid not in self.spectrum_uids_meeting_selection_criteria:
                continue
            comet_outputs = CometOutputs.standardized_comet_outputs(
                out_dir=self.hybrid_run_scan_results_dir,
                decoy_search=0,
                file_root=psm.sample,
                scan_min=psm.scan,
                scan_max=psm.scan,
            )
            expected_outputs.add(comet_outputs.target)
        return expected_outputs

    def run_native_assign_confidence(self, overwrite: bool = False):
        logger.info(f"Running `crux assign-confidence` for config {self.name}")
        native_target_txts = []
        for output in self.expected_native_comet_outputs:
            assert (
                output.target.exists()
            ), f"Expected target txt does not exist: {output.target}"
            native_target_txts.append(output.target)
        Crux().run_assign_confidence(
            target_txts=native_target_txts,
            out_path=self.native_assign_confidence_path,
            overwrite=overwrite,
        )

    def create_spectra_plots(self):
        df, fig, axs = plot_spectra_histograms(
            spectra=self.spectra,
            # add_cnts=True
        )
        save_fig(
            path=self.name_dir / "spectra_histograms.png", title=self.name, fig=fig
        )

    def create_native_run_plots(self):
        fig, axs = fig_setup(nrows=3, ncols=1)
        self.native_comet_run.create_top_target_vs_decoy_scatterplot(ax=axs[0])
        self.native_comet_run.create_top_target_and_decoy_xcorr_plot(ax=axs[1])
        self.native_comet_run.create_targets_xcorr_range_plot(ax=axs[2])
        save_fig(fig=fig, path=self.native_run_dir / "xcorr_plots.png", title=self.name)

        # Protein abundance plot
        fig, axs = fig_setup(h=10, w=8)
        self.protein_abundance_plot(
            ax=axs[0], q_threshold=DEFAULT_Q_THRESHOLD, top_n_prots_to_show=30
        )
        save_fig(
            fig=fig,
            path=self.native_run_dir / "protein_abundance_plot.png",
            title=self.name,
        )

        # Pair plot
        p = (
            self.native_comet_run.create_xcorr_target_decoy_diff_vs_peptide_len_jointplot()
        )
        p.fig.suptitle(self.name)
        p.savefig(self.native_run_dir / "xcorr_top_target_vs_top_decoy.png")

    def to_dict(self):
        d = asdict(self, dict_factory=path_aware_dict_factory)
        d["hybrid_run_params"] = d["hybrid_run_params"].model_dump(mode="json")
        return d

    def save(
        self,
        path: Union[Path, str],
    ):
        """Save config as a JSON"""
        save_dict(data=self.to_dict(), path=path)

    def run_param_medic(
        self,
        out_dir: Optional[str | Path] = None,
        crux_path: Optional[str] = MAC_CRUX_EXECUTABLE,
    ):
        if out_dir is None:
            out_dir = self.parent_output_dir / "param_medic"
            out_dir.mkdir(exist_ok=True, parents=True)
        for mzml_path in self.mzml_to_scans.keys():
            Mzml(path=mzml_path).run_param_medic(out_dir=out_dir, crux_path=crux_path)

    def run_native_comet_on_spectrum(
        self,
        mzml: Union[Path, str],
        scan: int,
        out_dir: Union[Path, str],
        decoy_search: Literal[0, 1, 2] = 0,
        crux_path: Optional[str | Path] = None,
    ) -> CometRun:
        comet_run = CometRun(
            mzml=mzml,
            fasta=self.hybrid_run_params.fasta,
            crux_comet_params=self.hybrid_run_params.crux_comet_params,
            out_dir=out_dir,
            decoy_search=decoy_search,
            scan_min=scan,
            scan_max=scan,
            num_threads=1,
        )
        comet_run.run_comet_and_keep_only_results(crux_path=crux_path)
        return comet_run

    def native_comet_run_on_all_spectra(
        self,
        crux_path: Optional[str | Path] = None,
        dry_run: bool = False,
    ) -> List[CometOutputs]:
        expected_outputs = []
        mzmls = list(self.mzml_to_scans.keys())
        for idx, (mzml, scans) in enumerate(self.mzml_to_scans.items()):
            if scans == ALL:
                comet_run = CometRun(
                    fasta=self.hybrid_run_params.fasta,
                    mzml=mzml,
                    crux_comet_params=self.hybrid_run_params.crux_comet_params,
                    decoy_search=2,
                    out_dir=self.native_run_dir,
                    file_root=Mzml.get_mzml_name(mzml=mzml),
                    dry_run=dry_run,
                )
                expected_outputs.append(comet_run.standardized_comet_outputs)
                if not dry_run:
                    logger.info(
                        f"Native Comet run on MZML {mzml.name} (MZML {idx+1} of {len(mzmls)})"
                    )
                    if comet_run.standardized_comet_outputs.target.exists():
                        logger.info(
                            "Looks like target TXT already exists so skipping..."
                        )
                        continue
                    comet_run.run_comet_and_keep_only_results(crux_path=crux_path)
            else:
                for scan in scans:
                    comet_run = CometRun(
                        fasta=self.hybrid_run_params.fasta,
                        mzml=mzml,
                        scan_min=scan,
                        scan_max=scan,
                        crux_comet_params=self.hybrid_run_params.crux_comet_params,
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
                            logger.info(
                                "Looks like target TXT already exists so skipping..."
                            )
                            continue
                        comet_run.run_comet_and_keep_only_results(
                            crux_path=crux_path,
                        )

        return expected_outputs

    def get_protein_abundance(
        self,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
    ) -> ProteinAbundance:
        return ProteinAbundance.from_comet_psms(
            psms=self.native_assign_confidence_psms, q_threshold=q_threshold
        )

    def protein_abundance_plot(
        self,
        ax: Axes | None = None,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
        top_n_prots_to_show: Optional[int] = None,
    ):
        prot_ab = ProteinAbundance.from_comet_psms(
            psms=self.native_assign_confidence_psms, q_threshold=q_threshold
        )
        prot_ab.plot_sorted_prot_cnts(
            top_n_prots=top_n_prots_to_show, ax=ax, q_threshold=q_threshold
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
        fasta = Fasta(path=self.hybrid_run_params.fasta)
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

    def combine_hybrid_run_scan_outputs(self, overwrite: bool = True):
        logger.info(f"Combining Comet scan results for config: {self.name}")
        mzml_to_txts = defaultdict(list)
        for txt in self.expected_hybrid_run_spectrum_target_txts:
            assert txt.exists(), f"Expected txt does not exist: {txt}"
            mzml, _, _ = CometOutputs.parse_standardized_comet_txt(comet_txt=txt)
            mzml_to_txts[mzml].append(txt)

        logger.info(
            "It seems like all the expected hybrid results are present. Combining them now..."
        )
        for mzml, txts in mzml_to_txts.items():
            logger.info(f"Combining Comet scan results for MZML {mzml}...")
            out_path = Path(self.hybrid_run_dir) / self.get_combined_comet_output_name(
                mzml_name=mzml, psm_type=TARGET
            )
            if out_path.exists() and overwrite:
                logger.info(
                    f"Combined Comet output {out_path} already exists. Overwrite is True so overwriting..."
                )
                out_path.unlink()
            _ = Crux.combine_crux_comet_files(
                files=txts,
                out_path=out_path,
            )
        logger.info("Finished combining Comet scan results")

    def check_for_missing_scans(
        self, print_missing: bool = False, raise_error: bool = False
    ):
        missing_spectra_paths = list(self.missing_hybrid_run_scan_target_txts)
        logger.info(
            f"There are {len(missing_spectra_paths)} spectra missing hybrid Comet outputs. "
        )
        if print_missing and len(missing_spectra_paths) > 0:
            logger.info("Spectra with missing HS outputs:")
            logger.info("\n".join([str(p) for p in missing_spectra_paths]))
        if raise_error:
            raise RuntimeError("There are spectra with no HypedSearch outputs!")

    def get_spectra_with_no_hybrid_results(self) -> List[Spectrum]:
        missing_spectra_paths = list(self.missing_hybrid_run_scan_target_txts)
        logger.info(
            f"There are {len(missing_spectra_paths)} spectra missing hybrid Comet outputs."
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
                self.spectrum_uid_to_spectrum[
                    Spectrum.get_uid(sample=name, scan=int(start_str))
                ]
            )
        logger.info("Done gathering missing spectra.")
        return missing_spectra

    def hybrid_run_on_spectrum(
        self,
        spectrum: Spectrum,
        fasta_dir: Path,
        crux_path: Optional[str | Path] = None,
        delete_hybrids_fasta: bool = True,
    ):
        return hybrid_run_on_spectrum(
            spectrum=spectrum,
            params=self.hybrid_run_params,
            fasta_dir=fasta_dir,
            crux_path=crux_path,
            delete_hybrids_fasta=delete_hybrids_fasta,
            out_dir=self.hybrid_run_scan_results_dir,
        )

    # Static methods
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


@log_time(level=logging.INFO)
def hybrid_run_on_spectrum(
    spectrum: Spectrum,
    params: HybridRunParams,
    fasta_dir: Path,
    out_dir: Path,
    crux_path: Optional[str | Path] = None,
    delete_hybrids_fasta: bool = True,
    overwrite: bool = False,
) -> Tuple[Optional[CmdLineResult], CometRun, Optional[Dict]]:
    # Get params and constants
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    hybrids_fasta = fasta_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    comet_run = CometRun(
        fasta=hybrids_fasta,
        mzml=spectrum.mzml,
        crux_comet_params=params.crux_comet_params,
        out_dir=out_dir,
        decoy_search=2 if params.hybrid_decoy_competition else 0,
        scan_min=spectrum.scan,
        scan_max=spectrum.scan,
        num_threads=1,
    )
    if comet_run.standardized_comet_outputs.target.exists() and not overwrite:
        logger.info(
            f"Comet output for spectrum {spectrum.uid} already exists at {comet_run.standardized_comet_outputs.target} and overwrite is False so skipping..."
        )
        return None, comet_run, None

    # Form hybrids
    hybrid_seq_to_position_strs = params.form_hybrids(spectrum=spectrum)
    if len(hybrid_seq_to_position_strs) == 0:
        logger.debug(
            f"No hybrids found for spectrum {spectrum.uid}. Skipping but creating empty Comet outputs..."
        )
        outputs = comet_run.standardized_comet_outputs
        outputs.target.touch()
        if outputs.decoy:
            outputs.decoy.touch()
        return None, comet_run, None

    # Create FASTA containing the hybrid peptides
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    hybrids_fasta = fasta_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    if params.hybrid_decoy_competition:
        create_hybrids_fasta(
            hybrid_seqs=set(hybrid_seq_to_position_strs.keys()),
            output_fasta_path=hybrids_fasta,
            fasta_to_include=params.fasta,
        )
    else:
        create_hybrids_fasta(
            hybrid_seqs=set(hybrid_seq_to_position_strs.keys()),
            output_fasta_path=hybrids_fasta,
        )
    process = comet_run.run_comet_and_keep_only_results(crux_path=crux_path)

    if not check_if_file_is_empty(comet_run.standardized_comet_outputs.target):
        # Update the "protein" column of the Comet output to include the positions that the hybrid sequence appears
        logger.debug(
            "Updating 'protein' column of Comet output to include hybrid position info..."
        )
        hybrid_psm_df = CometPSM.from_txt(
            txt=comet_run.standardized_comet_outputs.target, as_df=True
        )
        hybrid_psm_df[PROTEIN] = hybrid_psm_df[PLAIN_PEPTIDE].apply(
            lambda seq: COMET_PROTEIN_SEPARATOR.join(hybrid_seq_to_position_strs[seq])
        )
        # Overwrite results file
        hybrid_psm_df.to_csv(
            comet_run.standardized_comet_outputs.target, sep="\t", index=False
        )

        # Delete hybrids FASTA to save space if desired
        if delete_hybrids_fasta:
            logger.debug("Deleting hybrids FASTA")
            os.remove(hybrids_fasta)

    return process, comet_run, hybrid_seq_to_position_strs


def create_hybrids_fasta(
    hybrid_seqs: Set[str],
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
        [Peptide(seq=seq, name=f"hybrid{idx}") for idx, seq in enumerate(hybrid_seqs)]
    )
    Fasta.write_fasta(peptides=prots, path=output_fasta_path)
    return prots


@log_params
def run_hypedsearch(
    config: Path,
    n_cores: int,
    crux_path: Optional[str | Path] = None,
    stop_on_fail: bool = False,
    run_in_parallel: bool = True,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)

    # Get spectra to run Hypedsearch on
    # missing_spectra = hs_config.get_spectra_with_no_hybrid_results()
    spectra = [
        sp
        for sp in hs_config.spectra
        if sp.precursor_charge <= hs_config.hybrid_run_params.max_precursor_charge
    ]

    # Run Hypedsearch
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create FASTA containing hybrids and run Comet
        logger.info(f"Running Hypedsearch with FASTA dir: {tmp_dir}")
        if run_in_parallel:
            logger.info("Running in parallel")
            process_partial = partial(
                hybrid_run_on_spectrum,
                params=hs_config.hybrid_run_params,
                fasta_dir=Path(tmp_dir),
                out_dir=hs_config.hybrid_run_scan_results_dir,
                crux_path=crux_path,
                delete_hybrids_fasta=True,
            )
            with ProcessPoolExecutor(max_workers=n_cores) as ex:
                future_to_spectrum = {
                    ex.submit(process_partial, spectrum): spectrum.uid
                    for spectrum in spectra
                }
                # Process results as they complete, catching failures
                results = []
                failed = []
                for _, future in enumerate(as_completed(future_to_spectrum)):
                    spectrum_uid = future_to_spectrum[future]
                    try:
                        results.append(future.result())
                    except Exception as e:
                        msg = f"Task failed for spectrum {spectrum_uid}: {e}"
                        failed.append((spectrum_uid, str(e)))
                        logger.info(msg)
                        if stop_on_fail:
                            # raise
                            sys.exit(msg)

            if len(failed) > 0:
                raise RuntimeError(
                    f"{len(failed)}/{len(future_to_spectrum)} tasks failed."
                )

        else:
            logger.info("Running serially")
            for spectrum in spectra:
                hybrid_run_on_spectrum(
                    spectrum=spectrum,
                    params=hs_config.hybrid_run_params,
                    fasta_dir=Path(tmp_dir),
                    crux_path=crux_path,
                    out_dir=hs_config.hybrid_run_scan_results_dir,
                )
    logger.info(f"Finished running HypedSearch on config {config.name}")
