import re
from collections import defaultdict
from dataclasses import Field, dataclass
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union
from venv import logger

import pandas as pd
import yaml
from matplotlib import pyplot as plt
from pydantic import BaseModel, model_validator
from typing_extensions import Self

from src.comet_utils import CometPSM, CometPSMs
from src.constants import (
    COMET_DIR,
    CRUX_PATH_IN_SINGULARITY,
    DECOY,
    DEFAULT_CRUX_PARAMS,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_CLUSTER_LENGTH,
    DEFAULT_MIN_CLUSTER_SUPPORT,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_Q_VALUE_THRESH,
    HS_PREFIX,
    HYBRID,
    NATIVE,
    RUN_COMET_SMK,
    RUN_HYPEDSEARCH_SMK,
    TARGET,
)
from src.crux import CometConfig, CometOutputs, Crux, get_expected_comet_outputs
from src.hybrids_via_clusters import form_spectrum_hybrids_via_clustering
from src.kmer_database import KmerDatabase, create_kmer_database
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
from src.protein_abundance import (
    get_most_common_proteins,
    get_protein_counts_from_comet_results,
)
from src.utils import (
    copy_file,
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


class HybridRunConfig(BaseModel):
    mzml_to_scans: Dict[Path, List[int]]
    out_dir: Path
    fasta: Path
    kmer_db: Path
    kmer_to_protein_map: Path
    crux_path: Path = COMET_DIR / "crux-4.3.Darwin.x86_64/bin/crux"
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL
    crux_comet_params: Path = DEFAULT_CRUX_PARAMS
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT
    log_dir: Path = Path("")
    num_peaks: int = 0

    def save(self, path: Union[str, Path]) -> None:
        data = self.model_dump(mode="json")  # Converts Paths to strings
        save_dict(data=data, path=path)

    @classmethod
    def from_yaml(cls, path: Path) -> "HybridRunConfig":
        return cls(**load_yaml(path))

    @property
    def expected_outputs(self):
        return get_expected_comet_outputs(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=Path(self.out_dir),
            decoy_search=2,
            psm_type=TARGET,
        )

    def create_config_for_fiji_node(
        self, node_data_dir: Path, out_path: Optional[Path] = None
    ) -> "HybridRunConfig":
        logger.info(
            f"Preparing files for Hypedsearch on Fiji node in directory {node_data_dir}..."
        )
        # Handle MZMLs
        new_mzml_to_scans = {}
        for mzml, scans in self.mzml_to_scans.items():
            new_mzml_path = str(node_data_dir / Path(mzml).name)
            new_mzml_to_scans[new_mzml_path] = scans
            copy_file(src=mzml, dest=new_mzml_path)

        # Handle files that need to be copied
        fiji_config = {}
        for attr in ["kmer_db", "fasta", "crux_comet_params", "kmer_to_protein_map"]:
            old_path = Path(getattr(self, attr))
            new_out_path = str(node_data_dir / old_path.name)
            fiji_config[attr] = new_out_path
            copy_file(src=old_path, dest=new_out_path)

        # Create the output directory
        new_out_dir = node_data_dir / Path(self.out_dir).name
        new_out_dir.mkdir(parents=True, exist_ok=True)

        # Create the new config
        fiji_node_config = self.__class__(
            mzml_to_scans=new_mzml_to_scans,
            out_dir=str(new_out_dir),
            fasta=fiji_config["fasta"],
            kmer_db=fiji_config["kmer_db"],
            kmer_to_protein_map=fiji_config["kmer_to_protein_map"],
            crux_path=CRUX_PATH_IN_SINGULARITY,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            crux_comet_params=fiji_config["crux_comet_params"],
            min_cluster_len=self.min_cluster_len,
            min_cluster_support=self.min_cluster_support,
            log_dir=self.log_dir,
        )
        if out_path is not None:
            logger.info(f"Saving config to {out_path}...")
            fiji_node_config.save(path=out_path)
        return fiji_node_config


class HybridPSM(BaseModel):
    psm: CometPSM
    hybrid: HybridPeptide


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
    log_dir: Path = Path("")
    num_peaks: int = 0

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
    def kmer_to_protein_map_path(self) -> Path:
        return self.name_dir / "kmer_to_proteins_map.json"

    @property
    def kmer_to_protein_map(self) -> Dict[str, List[str]]:
        return load_json(self.kmer_to_protein_map_path)

    @property
    def native_run_smk_config(self) -> Path:
        return self.native_run_dir / "native_run.smk.config.json"

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
        return self.hybrid_run_dir / "hybrid_run.smk.config.json"

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
        for mzml, scans in self.mzml_to_scans.items():
            if scans == [0]:
                logger.info(
                    f"Updating scans for {mzml} from [0] to all of its scans..."
                )
                self.mzml_to_scans[mzml] = Mzml(mzml=mzml).scans

    def native_comet_snakemake_config_and_cmd(self) -> CometConfig:
        # Create config for snakemake
        mzml_to_scans = {mzml: [0] for mzml in self.mzml_to_scans.keys()}
        config = CometConfig(
            mzml_to_scans=mzml_to_scans,
            crux_comet_params=self.crux_comet_params,
            decoy_search=2,
            fasta=self.fasta,
            out_dir=self.native_run_dir,
            crux_path=self.crux_path,
        )
        config.save(path=self.native_run_smk_config)

        # Print command to run snakemake
        cmd = config.run_comet_snakemake_command(config_path=self.native_run_smk_config)
        logger.info(
            f"Saved native Comet run snakemake config to {self.native_run_smk_config}. To run snakemake, use this command:\n{cmd}"
        )
        return (cmd, config)

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
            target_txts.append(comet_outputs.target)

    def native_assign_confidence(self):
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
            kmer_to_protein_map=self.kmer_to_protein_map_path,
            crux_path=self.crux_path,
            peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
            precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
            crux_comet_params=self.crux_comet_params,
            min_cluster_len=self.min_cluster_len,
            min_cluster_support=self.min_cluster_support,
            log_dir=self.log_dir,
            num_peaks=self.num_peaks,
        )
        hybrid_run_config.save(path=self.hybrid_run_smk_config)
        logger.info(
            f"Saved hybrid Comet run snakemake config to {self.hybrid_run_smk_config}. To run snakemake, use this command:\n"
            f"snakemake -s {RUN_HYPEDSEARCH_SMK} --configfile {self.hybrid_run_smk_config} ..."
        )

    def psms(self, psm_type: Literal[NATIVE, HYBRID]) -> CometPSMs:
        if psm_type == NATIVE:
            psms = CometPSMs(
                psms=CometPSM.from_txt(txt=self.native_assign_confidence_txt)
            )
        elif psm_type == HYBRID:
            psms = CometPSM.from_txt(txt=self.hybrid_assign_confidence_txt)
            # Remove those PSMs that aren't hybrids
            psms = [psm for psm in psms if psm.is_hybrid]
            psms = CometPSMs(psms=psms)
        else:
            raise ValueError(
                f"Invalid psm_type: {psm_type}. Should be either {NATIVE} or {HYBRID}."
            )
        return psms

    def high_confidence_psms(
        self,
        psm_type: Literal[NATIVE, HYBRID],
        q_value_threshold: float = DEFAULT_Q_VALUE_THRESH,
    ) -> CometPSMs:
        psms = self.psms(psm_type=psm_type)
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
        prot_counts = get_protein_counts_from_comet_results(psms=psms)
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
            kmer_to_proteins_path=self.kmer_to_protein_map_path,
            db_path=self.kmer_db,
            min_k=min_k,
            max_k=max_k,
        )

    def combine_scan_results_and_run_assign_confidence(self):
        output_regex = r"^(?P<mzml>(.+?))\.comet\.(?P<scan>\d+)-(?P=scan)\.(?P<psm_type>target|decoy)\.txt$"

        # Get Comet outputs for each MZML
        mzml_files = defaultdict(lambda: {TARGET: [], DECOY: []})
        for hybrid_txt in self.hybrid_scan_results_dir.glob("*.txt"):
            match = re.match(output_regex, hybrid_txt.name)
            mzml_name = match.groupdict()["mzml"]
            psm_type = match.groupdict()["psm_type"]
            mzml_files[mzml_name][psm_type].append(hybrid_txt)

        target_txts = []
        for mzml in self.mzml_to_scans.keys():
            mzml_name = Mzml.get_mzml_name(mzml=mzml)

            # Combine targets
            out_path = self.hybrid_run_dir / f"{mzml_name}.comet.{TARGET}.txt"
            target_txts.append(out_path)
            _ = Crux.combine_crux_comet_files(
                files=mzml_files[mzml_name][TARGET],
                out_path=out_path,
            )

            # Combine decoys
            out_path = self.hybrid_run_dir / f"{mzml_name}.comet.{DECOY}.txt"
            _ = Crux.combine_crux_comet_files(
                files=mzml_files[mzml_name][DECOY],
                out_path=out_path,
            )

        # Run assign-confidence
        _ = Crux(path=self.crux_path).run_assign_confidence(
            target_txts=target_txts, out_path=self.hybrid_assign_confidence_txt
        )

    def protein_counts(
        self, q_value_threshold: float = DEFAULT_Q_VALUE_THRESH
    ) -> Dict[str, int]:
        psms = self.high_confidence_psms(
            psm_type=NATIVE, q_value_threshold=q_value_threshold
        )
        return get_protein_counts_from_comet_results(psms=psms.psms)

    def protein_count_plot(self, q_value_threshold: float = DEFAULT_Q_VALUE_THRESH):
        protein_counts = self.protein_counts(q_value_threshold=q_value_threshold)

        # Plot
        items = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
        keys, values = zip(*items)

        _, axs = fig_setup(h=8)
        ax = axs[0]
        set_title_axes_labels(
            ax=ax, title="Protein counts", xlabel="Protein", ylabel="Count"
        )
        finalize(axs)
        ax.scatter(range(len(keys)), values)
        ax.set_xticks(range(len(keys)), keys, rotation=90, fontsize=8)

    def save_hybrids(
        self,
        min_side_len: int,
        q_value_threshold: float = DEFAULT_Q_VALUE_THRESH,
        out_path: Optional[Path] = None,
    ) -> pd.DataFrame:
        hybrid_psms = self.high_confidence_psms(
            psm_type=HYBRID, q_value_threshold=q_value_threshold
        ).psms
        data = []
        for psm in hybrid_psms:
            hybrids = find_possible_hybrids(
                seq=psm.seq,
                kmer_to_protein_map=self.kmer_to_protein_map,
                min_side_len=min_side_len,
            )
            # hybrids
            for hybrid in hybrids:
                data.append({**psm.to_dict(), **hybrid.model_dump()})
        df = pd.DataFrame(data)
        if out_path is not None:
            df.to_csv(out_path, index=False)
            logger.info(f"Saved hybrids to {out_path}")
        return df


@dataclass
class HypedsearchOutputs:
    target: Path
    decoy: Path


def create_and_score_hybrids_for_spectrum(
    spectrum: Spectrum,
    kmer_db: KmerDatabase,
    protein_name_to_seq_map: Dict[str, str],
    kmer_to_proteins_map: Dict[str, List[str]],
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
    min_cluster_len: int,
    min_cluster_support: int,
    crux_comet_params: Path,
    out_dir: Path,
    fasta: Union[str, Path],
    crux_path: Path,
    num_peaks: int = 0,
) -> CometOutputs:
    # Peak filtering
    if num_peaks > 0:
        logger.info(f"Filtering to top {num_peaks} peaks...")
        spectrum.filter_to_top_n_peaks(n=num_peaks)
    # Form hybrids
    logger.info(f"Forming hybrids...")
    seq_to_hybrids = form_spectrum_hybrids_via_clustering(
        spectrum=spectrum,
        kmer_db=kmer_db,
        protein_name_to_seq_map=protein_name_to_seq_map,
        kmer_to_proteins_map=kmer_to_proteins_map,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        min_cluster_len=min_cluster_len,
        min_cluster_support=min_cluster_support,
    )
    # Create FASTA containing hybrids and run Comet
    logger.info("Creating FASTA containing hybrids...")
    mzml_name = Mzml.get_mzml_name(mzml=spectrum.mzml)
    fasta_containing_hybrids_path = out_dir / f"{mzml_name}.{spectrum.scan}.fasta"
    create_hybrids_fasta(
        hybrid_seqs=seq_to_hybrids.keys(),
        output_fasta_path=fasta_containing_hybrids_path,
        fasta_to_include=fasta,
    )
    logger.info("Running Comet on hybrid-containing FASTA...")
    outputs = Crux(path=crux_path).run_comet(
        mzml=spectrum.mzml,
        fasta=fasta_containing_hybrids_path,
        crux_comet_params=crux_comet_params,
        decoy_search=2,  # always want to run this with decoy search on
        out_dir=out_dir,
        file_root=mzml_name,
        scan_min=spectrum.scan,
        scan_max=spectrum.scan,
    )
    logger.info("Comet run complete. Deleting hybrid-containing FASTA...")
    fasta_containing_hybrids_path.unlink()
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
    prot_counts = get_protein_counts_from_comet_results(psms=psms)
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
):
    out_dir.mkdir(parents=True, exist_ok=True)
    config = CometConfig(
        mzml_to_scans={mzml: [0] for mzml in mzmls},
        crux_comet_params=crux_comet_params,
        decoy_search=2,
        fasta=fasta,
        out_dir=out_dir,
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
    seq: str, kmer_to_protein_map: Dict[str, List[str]], min_side_len: int
) -> List[HybridPeptide]:
    possible_hybrids = []
    for breakpoint in range(min_side_len, len(seq) - min_side_len + 1):
        left = seq[:breakpoint]
        right = seq[breakpoint:]
        if (left in kmer_to_protein_map) and (right in kmer_to_protein_map):
            possible_hybrids.append(
                HybridPeptide(
                    left_seq=left,
                    right_seq=right,
                    left_proteins=kmer_to_protein_map[left],
                    right_proteins=kmer_to_protein_map[right],
                )
            )
    return possible_hybrids


def filter_hybrids_by_side_length(
    hybrid_psms: List[CometPSM], kmers: Set, min_side_len: int
) -> List[CometPSM]:
    filtered_psms = []
    for psm in hybrid_psms:
        possible_hybrids = find_possible_hybrids(
            seq=psm.seq, kmer_to_protein_map=kmers, min_side_len=min_side_len
        )
        if len(possible_hybrids) > 0:
            filtered_psms.append(psm)
    return filtered_psms
