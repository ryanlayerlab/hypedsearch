import os
import re
import time
from collections import defaultdict
from dataclasses import Field, dataclass
from functools import cached_property
from pathlib import Path
from typing import Counter, Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import pandas as pd
import seaborn as sns
import yaml
from matplotlib import pyplot as plt
from pydantic import BaseModel, model_validator
from typing_extensions import Self

from src.blastp import run_blastp
from src.comet_utils import CometPSM
from src.constants import (
    COMET_DIR,
    CRUX_PATH_IN_SINGULARITY,
    DECOY,
    DEFAULT_CRUX_PARAMS,
    DEFAULT_MAX_ALLOWED_ION_CHARGE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_CLUSTER_LENGTH,
    DEFAULT_MIN_CLUSTER_SUPPORT,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_Q_VALUE_THRESH,
    GIT_REPO_DIR,
    HS_PREFIX,
    HYBRID,
    NATIVE,
    RUN_COMET_SMK,
    RUN_HYPEDSEARCH_SMK,
    TARGET,
)
from src.crux import CometConfig, CometOutputs, Crux, get_expected_comet_outputs
from src.dataclasses import CometPSMs, PSMScoring
from src.hybrids_via_clusters import form_spectrum_hybrids_via_clustering
from src.kmer_database import KmerDatabase, create_kmer_database
from src.mass_spectra import Mzml, Spectrum
from src.peptide_spectrum_comparison import PSM
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
from src.protein_abundance import (
    get_most_common_proteins,
    get_protein_counts_from_comet_psms,
)
from src.utils import (
    PathType,
    copy_file,
    flatten_list_of_lists,
    get_time_in_diff_units,
    load_json,
    load_yaml,
    save_dict,
    to_json,
    write_new_line_separated_file,
)


class HybridPeptide(BaseModel):
    left_seq: str
    right_seq: str
    left_proteins: List[str]
    right_proteins: List[str]

    @property
    def seq(self) -> str:
        return self.left_seq + self.right_seq


@dataclass
class HybridPSM:
    hybrids: List[HybridPeptide]
    psm: PSM
    comet_psm: CometPSM

    def to_dicts(self) -> List[Dict]:
        data = []
        for hybrid in self.hybrids:
            data.append(
                {
                    "mzml": self.comet_psm.sample,
                    "scan": self.comet_psm.scan,
                    "precursor_mz": self.psm.spectrum.precursor_mz,
                    "precursor_charge": self.psm.spectrum.precursor_charge,
                    "retention_time": self.psm.spectrum.retention_time,
                    "seq": hybrid.seq,
                    "left_seq": hybrid.left_seq,
                    "right_seq": hybrid.right_seq,
                    "left_proteins": hybrid.left_proteins,
                    "right_proteins": hybrid.right_proteins,
                    "xcorr": self.comet_psm.xcorr,
                    "q_value": self.comet_psm.q_value,
                    "ions_matched": self.comet_psm.ions_matched,
                    "ions_total": self.comet_psm.ions_total,
                    "prop_ions_matched": self.comet_psm.ions_matched
                    / self.comet_psm.ions_total,
                    "prefixes_supported": self.psm.prefixes_supported,
                    "prop_prefixes_supported": len(self.psm.prefixes_supported)
                    / len(hybrid.seq),
                    "suffixes_supported": self.psm.suffixes_supported,
                    "prop_suffixes_supported": len(self.psm.suffixes_supported)
                    / len(hybrid.seq),
                    "prop_intensity_supported": self.psm.prop_intensity_supported,
                }
            )
        return data


class HybridRunConfig(BaseModel):
    mzml_to_scans: Dict[Path, List[int]]
    out_dir: Path
    fasta: Path
    kmer_db: Path
    kmer_to_proteins_map: Path
    crux_path: Path = COMET_DIR / "crux-4.3.Darwin.x86_64/bin/crux"
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    crux_comet_params: Path = DEFAULT_CRUX_PARAMS
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE
    log_dir: Path = Path("")
    num_peaks: int = 0

    def save(self, path: Union[str, Path]) -> None:
        data = self.model_dump(mode="json")  # Converts Paths to strings
        save_dict(data=data, path=path)

    @classmethod
    def from_json(cls, path: Path) -> "HybridRunConfig":
        return cls(**load_json(path=path))

    @property
    def expected_outputs(self):
        return get_expected_comet_outputs(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=Path(self.out_dir),
            decoy_search=2,
            psm_type=TARGET,
        )


def combine_comet_scan_results(scan_results_dir: Path, out_dir: Path):
    output_regex = r"^(?P<mzml>(.+?))\.comet\.(?P<scan>\d+)-(?P=scan)\.(?P<psm_type>target|decoy)\.txt$"

    # Get Comet outputs for each MZML
    logger.info("Grouping Comet outputs by MZML...")
    mzml_names = set()
    mzml_files = defaultdict(lambda: {TARGET: [], DECOY: []})
    for hybrid_txt in scan_results_dir.glob("*.txt"):
        match = re.match(output_regex, hybrid_txt.name)
        mzml_name = match.groupdict()["mzml"]
        psm_type = match.groupdict()["psm_type"]
        mzml_files[mzml_name][psm_type].append(hybrid_txt)
        mzml_names.add(mzml_name)

    target_txts = []
    for mzml_name in mzml_names:
        # Combine targets
        logger.info(f"Combining Comet target outputs for {mzml_name}...")
        out_path = out_dir / f"{mzml_name}.comet.{TARGET}.txt"
        target_txts.append(out_path)
        _ = Crux.combine_crux_comet_files(
            files=mzml_files[mzml_name][TARGET],
            out_path=out_path,
        )

        # Combine decoys
        logger.info(f"Combining Comet decoy outputs for {mzml_name}...")
        out_path = out_dir / f"{mzml_name}.comet.{DECOY}.txt"
        _ = Crux.combine_crux_comet_files(
            files=mzml_files[mzml_name][DECOY],
            out_path=out_path,
        )


class HypedsearchConfig(BaseModel):
    name: str
    mzml_to_scans: Dict[Path, List[int]]
    out_dir: Path
    fasta: Path
    crux_path: Path = COMET_DIR / "crux-4.3.Darwin.x86_64/bin/crux"
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    crux_comet_params: Path = DEFAULT_CRUX_PARAMS
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    max_allowed_ion_charge: int = DEFAULT_MAX_ALLOWED_ION_CHARGE
    log_dir: Path = Path("")
    num_peaks: int = 0
    max_precursor_charge: Optional[int] = None

    @property
    def name_dir(self) -> Path:
        return self.out_dir / self.name

    @property
    def native_run_dir(self) -> Path:
        return self.name_dir / "native_run"

    @property
    def hybrid_run_dir(self) -> Path:
        return self.name_dir / "hybrid_run"

    @property
    def kmer_db(self) -> Path:
        return self.name_dir / "kmers.db"

    @property
    def kmer_to_proteins_map_path(self) -> Path:
        return self.name_dir / "kmer_to_proteins_map.json"

    @property
    def kmer_to_proteins_map(self) -> Dict[str, List[str]]:
        return load_json(self.kmer_to_proteins_map_path)

    @property
    def native_run_smk_config(self) -> Path:
        return self.native_run_dir / "native_run_config.json"

    @property
    def native_assign_confidence_txt(self) -> Path:
        return self.native_run_dir / "assign-confidence.txt"

    @property
    def hybrid_scan_results_dir(self) -> Path:
        return self.hybrid_run_dir / "scan_results"

    @property
    def hybrid_assign_confidence_txt(self) -> Path:
        return self.hybrid_run_dir / "assign-confidence.txt"

    @property
    def hybrid_run_smk_config(self) -> Path:
        return self.hybrid_run_dir / "hybrid_run_config.json"

    def missing_hybrid_outputs(self) -> List[str]:
        missing_target_outputs = get_expected_comet_outputs(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=Path(self.hybrid_scan_results_dir),
            decoy_search=2,
            psm_type=TARGET,
        )

        logger.info(
            f"Looking in {self.hybrid_scan_results_dir}... Missing results for {len(missing_target_outputs)} scans"
        )
        return missing_target_outputs

    def top_proteins_txt(self, top_n: int) -> Path:
        return self.name_dir / f"top_{top_n}_proteins.txt"

    @model_validator(mode="after")
    def post_init(self) -> Self:
        # Create required directories if they don't exist
        self.name_dir.mkdir(parents=True, exist_ok=True)
        self.native_run_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_scan_results_dir.mkdir(parents=True, exist_ok=True)

        # Set scans
        self.set_scans()
        return self

    def set_scans(self):
        if self.max_precursor_charge is None:
            for mzml, scans in self.mzml_to_scans.items():
                if scans == [0]:
                    logger.info(
                        f"Updating scans for {mzml} from [0] to all of its scans..."
                    )
                    self.mzml_to_scans[mzml] = Mzml(mzml=mzml).msn_scan_numbers
        else:
            for mzml, scans in self.mzml_to_scans.items():
                if scans == [0]:
                    logger.info(
                        f"Updating scans for {mzml} from [0] to all of its scans with precursor charge <= {self.max_precursor_charge}..."
                    )
                    spectra = list(
                        filter(
                            lambda spectrum: spectrum.precursor_charge
                            <= self.max_precursor_charge,
                            Mzml(mzml=mzml).ms2_spectra,
                        )
                    )
                    self.mzml_to_scans[mzml] = [s.scan for s in spectra]

    def to_json(self, path: Union[str, Path]) -> None:
        data = self.model_dump(mode="json")
        to_json(data=data, path=path)

    @classmethod
    def from_json(cls, path: Union[str, Path]) -> "HypedsearchConfig":
        return cls(**load_json(path))

    def create_native_run_comet_via_snakemake_config(
        self,
        out_path: Union[Path, str],
        num_threads_4_comet: int,
    ) -> CometConfig:
        # Create config for snakemake
        logger.debug(
            f"Creating and saving native Comet run snakemake config to {out_path}..."
        )
        mzml_to_scans = {mzml: [0] for mzml in self.mzml_to_scans.keys()}
        config = CometConfig(
            mzml_to_scans=mzml_to_scans,
            crux_comet_params=self.crux_comet_params,
            decoy_search=2,
            fasta=self.fasta,
            out_dir=self.native_run_dir,
            crux_path=self.crux_path,
            num_threads=num_threads_4_comet,
        )
        config.save(path=out_path)
        return config

    def create_native_run_comet_via_snakemake_config_and_print_run_cmd(
        self, num_threads_4_comet: int = 0
    ) -> CometConfig:
        comet_config = self.create_native_run_comet_via_snakemake_config(
            num_threads_4_comet=num_threads_4_comet, out_path=self.native_run_smk_config
        )
        print(
            f"Here's the Snakemake command for the native Comet run:\n"
            f"{comet_config.get_cmd_2_run_comet_via_snakemake(config_path=self.native_run_smk_config)}"
        )
        return comet_config

    def native_comet_run(self):
        # Run Comet
        crux = Crux(path=self.crux_path)
        for mzml in self.mzml_to_scans.keys():
            _ = crux.run_comet(
                mzml=mzml,
                fasta=self.fasta,
                crux_comet_params=self.crux_comet_params,
                decoy_search=2,
                out_dir=self.native_run_dir,
                file_root=Mzml.get_mzml_name(mzml=mzml),
            )

    def run_native_assign_confidence(self):
        # Assign confidence
        target_txts = [
            CometOutputs.standardized_comet_outputs(
                out_dir=self.native_run_dir,
                file_root=Mzml.get_mzml_name(mzml=mzml),
                decoy_search=2,
            ).target
            for mzml in self.mzml_to_scans.keys()
        ]
        crux = Crux(path=self.crux_path)
        crux.run_assign_confidence(
            target_txts=target_txts,
            out_path=self.native_assign_confidence_txt,
        )

    def create_hybrid_run_snakemake_config_and_cmd(self) -> HybridRunConfig:
        hybrid_run_config = HybridRunConfig(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=self.hybrid_scan_results_dir,
            fasta=self.fasta,
            kmer_db=self.kmer_db,
            kmer_to_proteins_map=self.kmer_to_proteins_map_path,
            crux_path=self.crux_path,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            crux_comet_params=self.crux_comet_params,
            min_cluster_len=self.min_cluster_len,
            min_cluster_support=self.min_cluster_support,
            log_dir=self.log_dir,
            num_peaks=self.num_peaks,
            max_allowed_ion_charge=self.max_allowed_ion_charge,
        )
        hybrid_run_config.save(path=self.hybrid_run_smk_config)
        logger.info(
            f"Saved hybrid Comet run snakemake config to {self.hybrid_run_smk_config}. To run snakemake, use this command:\n"
            f"snakemake -s {RUN_HYPEDSEARCH_SMK.relative_to(GIT_REPO_DIR)} --configfile {self.hybrid_run_smk_config} ..."
        )
        return hybrid_run_config

    def psms(self, psm_type: Literal[NATIVE, HYBRID]) -> CometPSMs:
        if psm_type == NATIVE:
            psms = CometPSMs(
                psms=CometPSM.from_txt(txt=self.native_assign_confidence_txt)
            )
        elif psm_type == HYBRID:
            psms = CometPSMs(
                psms=CometPSM.from_txt(txt=self.hybrid_assign_confidence_txt)
            )

            # Remove those PSMs that aren't hybrids
            psms = remove_non_hybrid_psms(psms=psms.psms)
            psms = CometPSMs(psms=psms)
        else:
            raise ValueError(
                f"Invalid psm_type: {psm_type}. Should be either {NATIVE} or {HYBRID}."
            )
        return psms

    def high_confidence_native_psms(self, q_value_threshold: float) -> CometPSMs:
        return self.native_psms.get_high_confidence_psms(
            q_value_threshold=q_value_threshold
        )

    def high_confidence_hybrid_psms(
        self, q_value_threshold: float, remove_high_confidence_natives: bool = True
    ) -> CometPSMs:
        psms = self.hybrid_psms
        if remove_high_confidence_natives:
            logger.info(
                "Removing hybrid PSMs for which there's a high-confidence native PSM..."
            )
            high_conf_native_psms = self.high_confidence_native_psms(
                q_value_threshold=q_value_threshold
            )
            native_scans = set(
                (psm.sample, psm.scan) for psm in high_conf_native_psms.psms
            )
            psms = CometPSMs(
                psms=[
                    psm
                    for psm in psms.psms
                    if (psm.sample, psm.scan) not in native_scans
                ]
            )
        else:
            logger.info(
                "NOT removing hybrid PSMs for which there's a high-confidence native PSM..."
            )
        return psms.get_high_confidence_psms(q_value_threshold=q_value_threshold)

    def create_kmer_database_from_top_n_proteins(
        self,
        top_n_proteins: int,
        q_value_threshold: float = DEFAULT_Q_VALUE_THRESH,
        min_k: int = DEFAULT_MIN_KMER_LEN,
        max_k: int = DEFAULT_MAX_KMER_LEN,
    ):
        # Get PSMs and filter to those that pass the q-value threshold
        psms = self.high_confidence_psms(
            q_value_threshold=q_value_threshold, psm_type=NATIVE
        ).psms

        # Get most abundant proteins
        prot_counts = get_protein_counts_from_comet_psms(psms=psms)
        most_common_proteins = get_most_common_proteins(
            protein_counts=prot_counts, top_n=top_n_proteins
        )
        logger.info(f"Top {top_n_proteins} proteins:\n{most_common_proteins}")
        _ = write_new_line_separated_file(
            lines=most_common_proteins, path=self.top_proteins_txt(top_n=top_n_proteins)
        )
        _ = create_kmer_database(
            proteins=Fasta(path=self.fasta).get_proteins_by_name(
                protein_names=most_common_proteins
            ),
            kmer_to_proteins_path=self.kmer_to_proteins_map_path,
            db_path=self.kmer_db,
            min_k=min_k,
            max_k=max_k,
        )

    def combine_hybrid_comet_scan_results(self):
        combine_comet_scan_results(
            scan_results_dir=self.hybrid_scan_results_dir, out_dir=self.hybrid_run_dir
        )

    def hybrid_assign_confidence(self):
        # Get target txts
        target_txts = [
            self.hybrid_run_dir / f"{Mzml.get_mzml_name(mzml=mzml)}.comet.{TARGET}.txt"
            for mzml in self.mzml_to_scans.keys()
        ]
        # Run assign-confidence
        _ = Crux(path=self.crux_path).run_assign_confidence(
            target_txts=target_txts, out_path=self.hybrid_assign_confidence_txt
        )

    def protein_counts(
        self, q_value_threshold: float = DEFAULT_Q_VALUE_THRESH
    ) -> Dict[str, int]:
        psms = self.native_psms.get_high_confidence_psms(
            q_value_threshold=q_value_threshold
        )
        return get_protein_counts_from_comet_psms(psms=psms.psms)

    def protein_count_plot(
        self,
        ax: plt.Axes,
        q_value_threshold: float = DEFAULT_Q_VALUE_THRESH,
        top_n_most_common: Optional[int] = None,
    ):
        protein_counts = self.protein_counts(q_value_threshold=q_value_threshold)

        # Plot
        items = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
        if top_n_most_common is not None:
            items = items[:top_n_most_common]
        keys, values = zip(*items)

        ax.scatter(range(len(keys)), values)
        ax.set_xticks(range(len(keys)), keys, rotation=90, fontsize=8)
        set_title_axes_labels(
            ax=ax,
            # title="Protein counts",
            xlabel="Protein",
            ylabel="PSM counts",
        )

    def get_hybrids(
        self,
        min_side_len: int,
        q_value_threshold: float,
    ) -> pd.DataFrame:
        psms = self.psms(psm_type="hybrid").psms
        psms = [psm for psm in psms if psm.q_value <= q_value_threshold]
        return get_hybrids_for_comet_psms(
            psms=psms,
            min_side_len=min_side_len,
            kmer_to_proteins_map=self.kmer_to_proteins_map,
        )

    def decoy_psms(self, psm_type: Literal[NATIVE, HYBRID]) -> CometPSMs:
        if psm_type == NATIVE:
            folder = self.native_run_dir
        elif psm_type == HYBRID:
            folder = self.hybrid_run_dir
        else:
            raise ValueError(
                f"Invalid psm_type: {psm_type}. Should be either {NATIVE} or {HYBRID}."
            )
        psms = flatten_list_of_lists(
            [CometPSM.from_txt(txt=txt) for txt in folder.glob(f"*{DECOY}.txt")]
        )
        psms = list(filter(lambda psm: psm.num == 1, psms))
        return CometPSMs(psms=psms)

    @cached_property
    def native_decoy_psms(self) -> CometPSMs:
        return self.decoy_psms(psm_type=NATIVE)

    @cached_property
    def hybrid_decoy_psms(self) -> CometPSMs:
        return self.decoy_psms(psm_type=HYBRID)

    @cached_property
    def native_psms(self) -> CometPSMs:
        return self.psms(psm_type=NATIVE)

    @cached_property
    def hybrid_psms(self) -> CometPSMs:
        return self.psms(psm_type=HYBRID)


@dataclass
class HypedsearchOutputs:
    target: Path
    decoy: Path


def remove_non_hybrid_psms(psms: List[CometPSM]) -> List[CometPSM]:
    logger.debug("Remove non-hybrid PSMs...")
    return [psm for psm in psms if psm.is_hybrid]


def get_hybrids_for_comet_psms(
    psms: List[CometPSM],
    min_side_len: int,
    kmer_to_proteins_map: Dict[str, List[str]],
) -> List[HybridPSM]:
    data = []
    for psm in psms:
        hybrids = find_possible_hybrids(
            seq=psm.seq,
            kmer_to_proteins_map=kmer_to_proteins_map,
            min_side_len=min_side_len,
        )
        if len(hybrids) > 0:
            data.append(HybridPSM(hybrids=hybrids, psm=psm))
    return data


def create_and_score_hybrids_for_spectrum(
    spectrum: Spectrum,
    kmer_db: KmerDatabase,
    protein_name_to_seq_map: Dict[str, str],
    kmer_to_proteins_map: Dict[str, List[str]],
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
    min_cluster_len: int,
    min_cluster_support: int,
    max_allowed_ion_charge: int,
    crux_comet_params: Path,
    out_dir: Path,
    fasta: Union[str, Path],
    crux_path: Path,
    num_peaks: int = 0,
) -> CometOutputs:
    start_time = time.perf_counter()
    # Peak filtering
    if num_peaks > 0:
        logger.info(f"Filtering to top {num_peaks} peaks...")
        spectrum.filter_to_top_n_peaks(n=num_peaks)
    # Form hybrids
    seq_to_hybrids = form_spectrum_hybrids_via_clustering(
        spectrum=spectrum,
        kmer_db=kmer_db,
        protein_name_to_seq_map=protein_name_to_seq_map,
        kmer_to_proteins_map=kmer_to_proteins_map,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        min_cluster_len=min_cluster_len,
        min_cluster_support=min_cluster_support,
        max_allowed_ion_charge=max_allowed_ion_charge,
    )
    # Create FASTA containing hybrids and run Comet
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    fasta_containing_hybrids_path = out_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    create_hybrids_fasta(
        hybrid_seqs=seq_to_hybrids.keys(),
        output_fasta_path=fasta_containing_hybrids_path,
        fasta_to_include=fasta,
    )

    # Run Comet
    outputs = Crux(path=crux_path).run_comet(
        mzml=spectrum.mzml,
        fasta=fasta_containing_hybrids_path,
        crux_comet_params=crux_comet_params,
        decoy_search=2,  # always want to run this with decoy search on
        out_dir=out_dir,
        file_root=mzml_name,
        scan_min=spectrum.scan,
        scan_max=spectrum.scan,
        num_threads=1,
    )
    # Comet run complete. Deleting hybrid-containing FASTA...
    fasta_containing_hybrids_path.unlink()
    duration = time.perf_counter() - start_time
    logger.info(
        f"Running HS on spectrum ({spectrum.sample}, {spectrum.scan}) took {get_time_in_diff_units(duration)}"
    )
    return HypedsearchOutputs(target=outputs.target, decoy=outputs.decoy)


def create_kmer_database_from_top_n_proteins(
    psms: List[CometPSM],
    top_n_proteins: int,
    fasta: Path,
    db_path: Path,
    top_proteins_txt: Path,
    kmer_to_protein_map_path: Path,
    min_k: int,
    max_k: int,
):
    prot_counts = get_protein_counts_from_comet_psms(psms=psms)
    most_common_proteins = get_most_common_proteins(
        protein_counts=prot_counts, top_n=top_n_proteins
    )
    logger.info(f"Top {top_n_proteins} proteins:\n{most_common_proteins}")
    _ = write_new_line_separated_file(lines=most_common_proteins, path=top_proteins_txt)
    _ = create_kmer_database(
        kmer_to_proteins_path=kmer_to_protein_map_path,
        fasta=fasta,
        proteins=most_common_proteins,
        db_path=db_path,
        min_k=min_k,
        max_k=max_k,
    )


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

    Fasta.write_fasta(peptides=prots, out_path=output_fasta_path)
    return prots


def create_native_comet_run_config(
    out_dir: Path,
    mzmls: List[Path],
    crux_comet_params: Path,
    fasta: Path,
    crux_path: Path,
    out_path: Path,
    num_threads: int = 0,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    config = CometConfig(
        mzml_to_scans={mzml: [0] for mzml in mzmls},
        crux_comet_params=crux_comet_params,
        decoy_search=2,
        fasta=fasta,
        out_dir=out_dir,
        num_threads=num_threads,
    ).to_dict()
    config["crux_path"] = str(crux_path)
    to_json(
        data=config,
        path=out_path,
    )
    logger.info(f"Saved native Comet run config to {out_path}")
    logger.info(
        "To run Comet via snakemake, use the following command:\n"
        f"snakemake -s snakefiles/run_comet.smk --configfile {out_path} ..."
    )


def find_possible_hybrids(
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
class HSRunAnalysis:
    native_targets: List[CometPSM]
    native_assign_confidence: List[CometPSM]
    native_decoys: List[CometPSM]
    hybrid_targets: List[CometPSM]
    hybrid_assign_confidence: List[CometPSM]
    hybrid_decoys: List[CometPSM]
    kmer_to_proteins_map: Dict[str, List[str]]
    sample_scan_to_spectrum_map: Dict[Tuple[str, int], Spectrum]

    def __post_init__(self):
        self.hybrid_targets = remove_non_hybrid_psms(psms=self.hybrid_targets)
        self.hybrid_assign_confidence = remove_non_hybrid_psms(
            psms=self.hybrid_assign_confidence
        )
        assert all([psm.is_hybrid for psm in self.hybrid_targets])
        assert all([psm.is_hybrid for psm in self.hybrid_assign_confidence])
        assert all([psm.num == 1 for psm in self.native_targets])
        assert all([psm.num == 1 for psm in self.native_decoys])
        assert all([psm.num == 1 for psm in self.hybrid_targets])
        assert all([psm.num == 1 for psm in self.hybrid_decoys])
        assert all([psm.num == 1 for psm in self.native_assign_confidence])
        assert all([psm.num == 1 for psm in self.hybrid_assign_confidence])

    def protein_counts(self, native_q_value_threshold: float):
        return get_protein_counts_from_comet_psms(
            psms=list(
                filter(
                    lambda psm: psm.q_value <= native_q_value_threshold,
                    self.native_assign_confidence,
                )
            )
        )

    def xcorr_distributions(self):
        fig, axs = fig_setup(nrows=1, ncols=2)
        ax = axs[0]
        data = [psm.xcorr for psm in self.native_targets]
        _ = sns.histplot(
            data,
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            label=f"Native run targets (n={len(data)})",
            fill=False,
        )
        data = [psm.xcorr for psm in self.native_decoys]
        _ = sns.histplot(
            data,
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            label=f"Native run decoys (n={len(data)})",
            fill=False,
        )
        max_native_decoy_xcorr = max(psm.xcorr for psm in self.native_decoys)
        ax.axvline(
            x=max_native_decoy_xcorr,
            color="black",
            linestyle="--",
            label=f"Max native decoy xcorr (={max_native_decoy_xcorr:.2f})",
        )

        ax = axs[1]
        data = [psm.xcorr for psm in self.hybrid_targets]
        _ = sns.histplot(
            data,
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            label=f"Hybrid run targets (only hybrids) (n={len(data)})",
            fill=False,
        )
        data = [psm.xcorr for psm in self.hybrid_decoys]
        _ = sns.histplot(
            data,
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            label=f"Hybrid run decoys (n={len(data)})",
            fill=False,
        )
        max_hybrid_decoy_xcorr = max(psm.xcorr for psm in self.hybrid_decoys)
        ax.axvline(
            x=max_hybrid_decoy_xcorr,
            color="grey",
            linestyle="--",
            label=f"Max hybrid decoy xcorr (={max_hybrid_decoy_xcorr:.2f})",
        )

        for ax in axs:
            set_title_axes_labels(
                ax=ax,
                ylabel="Density",
                xlabel="xcorr",
            )
        finalize(axs)
        return fig, axs

    def q_value_counter(self, psm_type: Literal[NATIVE, HYBRID]):
        if psm_type == NATIVE:
            return Counter(psm.q_value for psm in self.native_assign_confidence)
        elif psm_type == HYBRID:
            return Counter(psm.q_value for psm in self.hybrid_assign_confidence)
        else:
            raise ValueError(f"Invalid psm_type: {psm_type}.")

    def xcorr_change(self, q_value_threshold: Optional[float] = None):
        # Get the spectra that have both a native and hybrid target PSM
        if q_value_threshold is None:
            native_spectra = set((psm.sample, psm.scan) for psm in self.native_targets)
            hybrid_spectra = set((psm.sample, psm.scan) for psm in self.hybrid_targets)
        else:
            native_spectra = set(
                (psm.sample, psm.scan)
                for psm in self.native_assign_confidence
                if psm.q_value <= q_value_threshold
            )
            hybrid_spectra = set(
                (psm.sample, psm.scan)
                for psm in self.hybrid_assign_confidence
                if psm.q_value <= q_value_threshold
            )
        shared_spectra = native_spectra.intersection(hybrid_spectra)

        # TODO: finish

    def set_decoy_scores(
        self, peak_to_ion_ppm_threshold: float, spectra_dir: Union[str, Path]
    ) -> List[PSMScoring]:
        logger.info("Getting 'PSM' objects for all decoys...")
        decoy_psms = PSM.from_comet_psms(
            comet_psms=self.native_decoys + self.hybrid_decoys,
            spectra_dir=spectra_dir,
            peak_to_ion_ppm_threshold=peak_to_ion_ppm_threshold,
        )

        logger.info("Calculating decoy scores and setting ECDFs...")
        ions_score = PSMScoring(name="prop_ions_matched")
        prefix_score = PSMScoring(name="prop_prefixes_supported")
        suffix_score = PSMScoring(name="prop_suffixes_supported")
        intensity_score = PSMScoring(name="prop_intensity_supported")
        for psm in decoy_psms:
            charge = psm.spectrum.precursor_charge
            ions_score.add_value(
                charge=charge,
                value=psm.prop_ions_matched,
            )
            prefix_score.add_value(charge=charge, value=psm.prop_prefixes_supported)
            suffix_score.add_value(charge=charge, value=psm.prop_suffixes_supported)
            intensity_score.add_value(charge=charge, value=psm.prop_intensity_supported)
        decoy_scores = [ions_score, prefix_score, suffix_score, intensity_score]
        for psm_score in decoy_scores:
            psm_score.set_ecdfs()
        self.decoy_scores = decoy_scores
        return decoy_scores

    @cached_property
    def max_decoy_xcorr(self) -> float:
        return max(psm.xcorr for psm in self.native_decoys + self.hybrid_decoys)

    def set_df(
        self,
        native_psm_q_value_threshold: float,
        min_side_len: int,
        peak_to_ion_ppm_tolerance: float,
        fasta: Union[str, Path],
    ) -> pd.DataFrame:
        # Remove hybrid PSMs that are also found as high-confidence native PSMs
        logger.info(
            f"Removing hybrid PSMs that have a high-confidence (q<={native_psm_q_value_threshold}) native PSM..."
        )
        native_psm_to_remove = set(
            (psm.sample, psm.scan)
            for psm in self.native_assign_confidence
            if psm.q_value <= native_psm_q_value_threshold
        )
        hybrid_comet_psms = self.hybrid_assign_confidence
        hybrid_comet_psms = [
            psm
            for psm in hybrid_comet_psms
            if (psm.sample, psm.scan) not in native_psm_to_remove
        ]

        # Get hybrid peptides for each hybrid PSM and create a HybridPSM object for each PSM
        logger.info(
            "For each hybrid PSM, get the PSMs possible explanatory hybrid sequences that might have produced it..."
        )
        hybrid_psms = []
        for comet_psm in hybrid_comet_psms:
            hybrids = find_possible_hybrids(
                seq=comet_psm.seq,
                kmer_to_proteins_map=self.kmer_to_proteins_map,
                min_side_len=min_side_len,
            )
            if len(hybrids) > 0:
                psm = PSM(
                    spectrum=self.sample_scan_to_spectrum_map[
                        (comet_psm.sample, comet_psm.scan)
                    ],
                    peptide=comet_psm.seq,
                    peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                )

                hybrid_psms.append(
                    HybridPSM(hybrids=hybrids, psm=psm, comet_psm=comet_psm)
                )

        df = pd.DataFrame(
            flatten_list_of_lists(list_of_lists=[psm.to_dicts() for psm in hybrid_psms])
        )

        # Add BLASTP results
        logger.info("Add BLASTP info...")
        blast_df = run_blastp(
            query_peptides=list(set(df["seq"].to_list())),
            fasta=fasta,
        )
        seq_to_blast_df = {}
        for seq, seq_df in blast_df.groupby("qseqid"):
            seq_to_blast_df[seq] = seq_df

        blast_coverages = []
        blast_percent_identities = []
        blast_num_mismatches = []
        for row_idx, row in df.iterrows():
            try:
                seq_df = seq_to_blast_df[row.seq]
                blast_row = seq_df.loc[seq_df["qcovhsp"].idxmax()]
                blast_coverages.append(blast_row["qcovhsp"])
                blast_percent_identities.append(blast_row["pident"])
                blast_num_mismatches.append(blast_row["mismatch"])
            except:
                blast_coverages.append(0)
                blast_percent_identities.append(0)
                blast_num_mismatches.append(len(row.seq))

        df["blast_coverage"] = blast_coverages
        df["blast_percent_identity"] = blast_percent_identities
        df["blast_num_mismatches"] = blast_num_mismatches

        self.df = df
        return df

    def add_protein_abundance_and_decoy_based_scores_to_df(
        self,
        native_psm_q_value_threshold: float,
        peak_to_ion_ppm_tolerance: float,
        spectra_dir: Union[str, Path],
    ):
        # Get ECDFs for each decoy score
        logger.info("Getting decoy-based score ECDFs...")
        self.set_decoy_scores(
            peak_to_ion_ppm_threshold=peak_to_ion_ppm_tolerance, spectra_dir=spectra_dir
        )

        # Add protein abundance and decoy-based scores to df
        protein_counts = self.protein_counts(
            native_q_value_threshold=native_psm_q_value_threshold
        )
        left_abs = []
        right_abs = []
        ion_probs = []
        prefix_probs = []
        suffix_probs = []
        intensity_probs = []

        ion_score = list(
            filter(lambda score: score.name == "prop_ions_matched", self.decoy_scores)
        )[0]
        prefix_score = list(
            filter(
                lambda score: score.name == "prop_prefixes_supported", self.decoy_scores
            )
        )[0]
        suffix_score = list(
            filter(
                lambda score: score.name == "prop_suffixes_supported", self.decoy_scores
            )
        )[0]
        intensity_score = list(
            filter(
                lambda score: score.name == "prop_intensity_supported",
                self.decoy_scores,
            )
        )[0]

        for row_idx, row in self.df.iterrows():
            # Protein abundances
            left_abs.append(max([protein_counts[prot] for prot in row.left_proteins]))
            right_abs.append(max([protein_counts[prot] for prot in row.right_proteins]))

            # Score probabilities (1-pvalue)
            ion_probs.append(
                float(
                    ion_score.ecdfs_by_charge[row.precursor_charge].cdf.evaluate(
                        row.prop_ions_matched
                    )
                )
            )
            prefix_probs.append(
                float(
                    prefix_score.ecdfs_by_charge[row.precursor_charge].cdf.evaluate(
                        row.prop_prefixes_supported
                    )
                )
            )
            suffix_probs.append(
                float(
                    suffix_score.ecdfs_by_charge[row.precursor_charge].cdf.evaluate(
                        row.prop_suffixes_supported
                    )
                )
            )
            intensity_probs.append(
                float(
                    intensity_score.ecdfs_by_charge[row.precursor_charge].cdf.evaluate(
                        row.prop_intensity_supported
                    )
                )
            )

        self.df["ion_prob"] = ion_probs
        self.df["prefix_prob"] = prefix_probs
        self.df["suffix_prob"] = suffix_probs
        self.df["intensity_prob"] = intensity_probs
        self.df["left_protein_ab_count"] = left_abs
        self.df["right_protein_ab_count"] = right_abs

        # Derived columns
        self.df["left_protein_ab"] = self.df["left_protein_ab_count"] / max(
            protein_counts.values()
        )
        self.df["right_protein_ab"] = self.df["right_protein_ab_count"] / max(
            protein_counts.values()
        )
        self.df["protein_ab_score"] = (
            self.df["left_protein_ab"] * self.df["right_protein_ab"]
        )
        self.df["min_prob"] = self.df.apply(
            lambda row: min(
                row["ion_prob"],
                row["prefix_prob"],
                row["suffix_prob"],
                row["intensity_prob"],
            ),
            axis=1,
        )


@click.command(
    name="combine-comet-scan-results",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=(
        "Given a directory containing the Comet output for each scan, combine them by mzML"
    ),
)
@click.option(
    "--scan_results_dir",
    "-srd",
    type=PathType(),
    required=True,
    help="Path to the directory containing scan results",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=True,
    help="Path to the output directory where the combined results for each mzML will be saved",
)
def cli_combine_comet_scan_results(scan_results_dir: Path, out_dir: Path):
    combine_comet_scan_results(scan_results_dir=scan_results_dir, out_dir=out_dir)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    cli.add_command(cli_combine_comet_scan_results)
    cli()
    cli()
    cli()
