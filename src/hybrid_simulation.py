import logging
import random
import shutil
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property, partial
from pathlib import Path
from typing import Counter, Dict, List, Literal, Optional, Tuple, Union

import click
import numpy as np
import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import NATIVE
from src.crux import CometRun, Crux
from src.hypedsearch import HybridRunParams, hybrid_run_on_spectrum
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.peptides_and_ions import Fasta, Fasta2MFMIndex, Peptide
from src.psm import CometPSM
from src.utils import PathType, load_json, log_params, move_file, setup_logger, to_json

logger = logging.getLogger(__name__)

DIR = Path(__file__).parent
HALF = "half"
RANDOM = "random"
LEFT_SEQ = "left_hybrid_seq"
RIGHT_SEQ = "right_hybrid_seq"
DEFAULT_HYBRID_SIM_CONFIG_NAME = "hybrid.simulation.config.json"
DENATIVIZED = "denativized"
# CUT_METHODS =
# CUT_METHODS_TYPE = Literal[HALF, RANDOM, ]


def cut_seq_into_hybrid(
    seq: str, min_side_len: int, method: Union[RANDOM, float, int] = 0.5
):
    assert len(seq) >= 2 * min_side_len
    if method == RANDOM:
        internal_seq = seq[min_side_len:-min_side_len]
        cut_idx = min_side_len + random.randint(0, len(internal_seq))
        left_hy_seq = seq[:cut_idx]
        right_hy_seq = seq[cut_idx:]
    elif isinstance(method, int):
        assert (
            0 < method < len(seq)
        ), "If method is an int, it should be between 0 and the length of the sequence."
        cut_idx = method
        left_hy_seq = seq[:cut_idx]
        right_hy_seq = seq[cut_idx:]
    elif isinstance(method, float):
        assert 0 < method < 1, "If method is a float, it should be between 0 and 1."
        cut_idx = int(len(seq) * method)
        left_hy_seq = seq[:cut_idx]
        right_hy_seq = seq[cut_idx:]
    else:
        raise ValueError(
            f"Unknown cut method: {method}. Allowed values are '{RANDOM}' and '{HALF}'."
        )
    return left_hy_seq, right_hy_seq


def create_new_kmer_db_and_fasta_for_simulation_with_hybrid_within_prot(
    existing_kmer_db: Union[str, Path],
    new_kmer_db: Union[str, Path],
    fasta: Union[str, Path, Fasta],
    psm_seq: str,
    left_seq: str,
    right_seq: str,
    new_fasta: Union[str, Path],
):
    """ """
    if isinstance(fasta, (str, Path)):
        fasta = Fasta(path=fasta)
    kmer_db = KmerDatabase(db_path=existing_kmer_db)
    kmer_db_prots = fasta.get_proteins_by_name(names=kmer_db.proteins)
    new_prots = []
    prot_containing_seq = None
    num_prots_containing_seq = 0
    for pep in kmer_db_prots:
        new_pep = pep.model_copy(deep=True)
        if psm_seq in pep.seq:
            prot_containing_seq = new_pep
            num_prots_containing_seq += 1
        else:
            new_prots.append(new_pep)
    assert num_prots_containing_seq == 1
    logger.info(
        f"Number of kmer database proteins containing seq {psm_seq}: {num_prots_containing_seq}"
    )
    # Update sequence-containing protein to have hybrid sequence
    split_seq = prot_containing_seq.seq.split(psm_seq)
    assert len(split_seq) == 2
    prot_containing_seq.seq = (
        split_seq[0][: int(len(split_seq[0]) / 2)]
        + left_seq
        + split_seq[0][int(len(split_seq[0]) / 2) :]
        + split_seq[1][: int(len(split_seq[1]) / 2)]
        + right_seq
        + split_seq[1][int(len(split_seq[1]) / 2) :]
    )
    new_prots.append(prot_containing_seq)
    KmerDatabase.create_db(db_path=new_kmer_db, proteins=new_prots, overwrite=True)

    new_prot_names = [prot.name for prot in new_prots]
    for prot in fasta.proteins:
        if prot.name not in new_prot_names:
            new_prots.append(prot)
    Fasta.write_fasta(peptides=new_prots, path=new_fasta)


def create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots(
    existing_kmer_db: KmerDatabase,
    fasta: Fasta,
    left_aa_seq: str,
    right_aa_seq: str,
    hybridized_kmer_db: Union[str, Path],
    hybridized_fasta: Union[str, Path],
):
    logger.info("Creating hybridized k-mer database and FASTA for simulation")
    # Constants
    aa_seq = left_aa_seq + right_aa_seq
    left_prot, right_prot = Peptide(seq=left_aa_seq, name=LEFT_SEQ), Peptide(
        seq=right_aa_seq, name=RIGHT_SEQ
    )
    # Validation and get the protein that contains the given sequence
    prots_containing_seq = fasta.proteins_that_contain_seqs(
        seqs=[left_aa_seq + right_aa_seq]
    )[aa_seq]
    assert (
        len(prots_containing_seq) == 1
    ), "Only one protein is expected to contain the hybrid sequence"
    prot_containing_seq_name = list(prots_containing_seq)[0]
    # Turn the protein into a Peptide object
    prot_containing_seq_with_seq_removed = Peptide(  # turn the prot
        name=f"{prot_containing_seq_name}|hybridized",
        seq=fasta.protein_name_to_seq_map[prot_containing_seq_name].replace(aa_seq, ""),
    )

    # Re-create the k-mer database
    kmer_db_prots = [left_prot, right_prot]
    if prot_containing_seq_name in existing_kmer_db.proteins:
        # Handle case when the k-mer DB contains the protein that contains the sequence
        prot_names = deepcopy(existing_kmer_db.proteins)
        prot_names.remove(prot_containing_seq_name)
        kmer_db_prots.extend(fasta.get_proteins_by_name(names=prot_names))
        kmer_db_prots.append(prot_containing_seq_with_seq_removed)
    else:
        kmer_db_prots.extend(
            fasta.get_proteins_by_name(names=existing_kmer_db.proteins)
        )

    # Create new k-mer database and FASTA
    KmerDatabase.create_db(
        db_path=hybridized_kmer_db,
        proteins=kmer_db_prots,
        overwrite=True,
        min_k=existing_kmer_db.min_k,
        max_k=existing_kmer_db.max_k,
    )
    new_prots_for_fasta = [
        prot for prot in fasta.proteins if prot.name != prot_containing_seq_name
    ] + [
        prot_containing_seq_with_seq_removed,
        left_prot,
        right_prot,
    ]
    Fasta.write_fasta(peptides=new_prots_for_fasta, path=hybridized_fasta)


def validate_hybrid_for_hybrid_simulation(
    left_seq: str,
    right_seq: str,
    fasta: Fasta,
    min_side_len: int,
) -> bool:
    """
    Make sure the given PSM is easy to turn from a native to a hybrid which means:
    1) The PSM sequence appears in only one protein
    2) The PSM sequence appears only once in that protein
    If the PSM satisifies these conditions, return True, else False
    """
    # Constants and setup
    aa_seq = left_seq + right_seq
    # Make sure the sequence appears in only one protein
    validation_failure_msg = (
        f"Validation failed for sequence {aa_seq} in FASTA {fasta.path}."
    )
    seq_containing_prots = list(fasta.proteins_that_contain_seqs([aa_seq])[aa_seq])
    if len(seq_containing_prots) != 1:
        logger.info(
            f"{validation_failure_msg} Sequence appears in multiple proteins: {seq_containing_prots}"
        )
        return False

    # Make sure sequence appears only once in that protein
    prot_name = seq_containing_prots[0]
    if fasta.protein_name_to_seq_map[prot_name].count(aa_seq) != 1:
        logger.info(
            f"{validation_failure_msg} Sequence appears in one protein but appears > 1 times in {prot_name}."
        )
        return False

    # Make sure sequence each side length is >= min_side_len
    if len(left_seq) < min_side_len or len(right_seq) < min_side_len:
        logger.info(
            f"{validation_failure_msg} One of the sides of the hybrid sequence is shorter than the minimum side length of {min_side_len}. Left seq: {left_seq}, right seq: {right_seq}."
        )
        return False
    return True


class HybridSimulationOutput(BaseModel):
    native_txt: Path
    denativized_txt: Path
    denativeized_hypedsearch_txt: Path


class HybridSimulationExperiment(BaseModel):
    hybrid_run_params: HybridRunParams
    scan_to_left_right_seq: Dict[int, Tuple[str, str]]
    mzml: Mzml
    parent_out_dir: Path

    @field_validator("parent_out_dir", mode="before")
    @classmethod
    def create_output_dir(cls, v):
        path = Path(v)
        path.mkdir(parents=True, exist_ok=True)
        return path

    @cached_property
    def native_out_dir(self) -> Path:
        d = self.parent_out_dir / "native_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @cached_property
    def hybrid_out_dir(self) -> Path:
        d = self.parent_out_dir / "denativized_hybrid_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def n(self):
        return len(self.scan_to_left_right_seq)

    @property
    def native_combined_txt(self) -> Path:
        return self.parent_out_dir / "combined_native_psms.txt"

    @property
    def hybrid_combined_txt(self) -> Path:
        return self.parent_out_dir / "combined_hybrid_psms.txt"

    @cached_property
    def native_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.native_combined_txt)

    @cached_property
    def hybrid_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.hybrid_combined_txt)

    @property
    def min_side_len_dict(self) -> Dict[int, int]:
        return dict(
            Counter(
                min([len(left_seq), len(right_seq)])
                for (left_seq, right_seq) in self.scan_to_left_right_seq.values()
            )
        )

    @classmethod
    def create_experiment(
        cls,
        hybrid_run_params: str | Path,
        mzml: str | Path,
        psms: List[CometPSM],
        parent_out_dir: str | Path,
        cut_method: Union[RANDOM, float, int] = 0.5,
        min_side_len: int = 3,
    ):
        assert set(Counter(psm.scan for psm in psms).values()) == set(
            [1]
        ), "Each scan should have exactly one PSM associated with it"
        hybrid_run_params = HybridRunParams.load(path=hybrid_run_params)
        fasta = Fasta(path=hybrid_run_params.fasta)
        scan_to_left_right_seqs = {}
        for psm in psms:
            left_seq, right_seq = cut_seq_into_hybrid(
                seq=psm.seq,
                method=cut_method,
                min_side_len=min_side_len,
            )
            if validate_hybrid_for_hybrid_simulation(
                left_seq=left_seq,
                right_seq=right_seq,
                fasta=fasta,
                min_side_len=min_side_len,
            ):
                scan_to_left_right_seqs[psm.scan] = (left_seq, right_seq)
            else:
                logger.info(
                    f"PSM (scan={psm.scan}) with sequence {psm.seq} failed validation. Ignoring this PSM"
                )
        parent_out_dir = Path(parent_out_dir)
        parent_out_dir.mkdir(parents=True, exist_ok=True)
        exp = cls(
            hybrid_run_params=hybrid_run_params,
            scan_to_left_right_seq=scan_to_left_right_seqs,
            mzml=Mzml(path=mzml),
            parent_out_dir=parent_out_dir,
        )

        # Save config to parent output directory
        exp.save(path=parent_out_dir / DEFAULT_HYBRID_SIM_CONFIG_NAME)
        return exp

    @classmethod
    def load(cls, path: str | Path) -> "HybridSimulationExperiment":
        data = load_json(path=path)
        if isinstance(data["hybrid_run_params"], (str, Path)):
            data["hybrid_run_params"] = HybridRunParams.load(
                path=data["hybrid_run_params"]
            )
        data["mzml"] = Mzml(path=data["mzml"])
        return cls(**data)

    def save(self, path: str | Path):
        data = self.model_dump(mode="json")
        data["mzml"] = str(self.mzml.path)
        to_json(
            data=data,
            path=path,
        )

    def run_hybrid_simulation_on_spectrum(
        self,
        scan: int,
        left_seq: str,
        right_seq: str,
        crux_path: Optional[str | Path] = None,
    ):
        run_hybrid_simulation_on_spectrum(
            scan=scan,
            left_seq=left_seq,
            right_seq=right_seq,
            mzml=self.mzml,
        )

    def run_experiment_in_parallel(
        self, n_cores: int, crux_path: Optional[str | Path] = None
    ):
        partial_fcn = partial(
            run_hybrid_simulation_on_spectrum,
            mzml=self.mzml,
            hybrid_run_params=self.hybrid_run_params,
            native_out_dir=self.native_out_dir,
            hybrid_out_dir=self.hybrid_out_dir,
            crux_path=crux_path,
        )
        with ProcessPoolExecutor(max_workers=n_cores) as ex:
            future_to_scan = {
                ex.submit(partial_fcn, scan, left_seq, right_seq): scan
                for scan, (left_seq, right_seq) in self.scan_to_left_right_seq.items()
            }
            for idx, future in enumerate(as_completed(future_to_scan)):
                try:
                    _ = future.result()
                    logger.info(
                        f"Finished scan {future_to_scan[future]} ({idx + 1} of {len(self.scan_to_left_right_seq)})"
                    )
                except Exception as e:
                    logger.warning(
                        f"Task failed for scan {future_to_scan[future]}: {e}"
                    )

    def combine_spectrum_txt_results(self):
        native_txts = list(self.native_out_dir.glob("*.txt"))
        Crux.combine_crux_comet_files(
            files=native_txts,
            out_path=self.native_combined_txt,
        )
        hybrid_txts = list(self.hybrid_out_dir.glob("*.txt"))
        Crux.combine_crux_comet_files(
            files=hybrid_txts, out_path=self.hybrid_combined_txt
        )

    def compare_native_psm_to_hybrid_psm(
        self, native_psm: CometPSM, hybrid_psm: CometPSM
    ) -> Dict:
        scan = native_psm.scan
        assert (
            hybrid_psm.scan == scan
        ), "PSMs should be from the same scan to be comparable"
        return {
            "scan": scan,
            "native_seq": native_psm.seq,
            "hybrid_seq": hybrid_psm.seq,
            "native_xcorr": native_psm.xcorr,
            "hybrid_xcorr": hybrid_psm.xcorr,
            "left_seq": self.scan_to_left_right_seq[scan][0],
            "right_seq": self.scan_to_left_right_seq[scan][1],
        }

    def compare_native_and_hybrid_psms(self) -> pd.DataFrame:
        # Get top native PSM for each scan
        native_scan_to_psms = defaultdict(list)
        for psm in self.native_psms:
            native_scan_to_psms[psm.scan].append(psm)
        native_scan_to_psms = {
            scan: sorted(psms, key=lambda psm: psm.num)
            for scan, psms in native_scan_to_psms.items()
        }
        hybrid_scan_to_psms = defaultdict(list)
        for psm in self.hybrid_psms:
            hybrid_scan_to_psms[psm.scan].append(psm)
        hybrid_scan_to_psms = {
            scan: sorted(psms, key=lambda psm: psm.num)
            for scan, psms in hybrid_scan_to_psms.items()
        }

        df = []
        for scan in native_scan_to_psms.keys():
            native_psm = native_scan_to_psms[scan][0]
            hybrid_psm = hybrid_scan_to_psms[scan][0]
            assert (native_psm.num == 1) and (hybrid_psm.num == 1)
            df.append(
                self.compare_native_psm_to_hybrid_psm(
                    native_psm=native_psm,
                    hybrid_psm=hybrid_psm,
                )
            )
        df = pd.DataFrame(df)
        df["equal"] = df["native_seq"] == df["hybrid_seq"]
        df["min_side_len"] = df.apply(
            lambda row: min([len(row.left_seq), len(row.right_seq)]), axis=1
        )
        return df


def run_hybrid_simulation_on_spectrum(
    scan: int,
    left_seq: str,
    right_seq: str,
    mzml: Mzml,
    hybrid_run_params: HybridRunParams,
    native_out_dir: Path,
    hybrid_out_dir: Path,
    crux_path: Optional[str | Path] = None,
):
    # Validation and load function constants
    assert (
        native_out_dir != hybrid_out_dir
    ), "Native and hybrid output directories must be different to avoid overwriting outputs since we want to maintain the MZML name in the output file name for both runs."
    fasta = Fasta(path=hybrid_run_params.fasta)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        hybridized_kmer_db = tmp_dir / "hybridized_kmers.db"
        hybridized_fasta = tmp_dir / "hybridized_proteins.fasta"
        hybridized_fasta_fm_index = tmp_dir / "hybridized_proteins.mfm"
        create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots(
            existing_kmer_db=hybrid_run_params.kmer_db,
            hybridized_kmer_db=hybridized_kmer_db,
            hybridized_fasta=hybridized_fasta,
            fasta=fasta,
            left_aa_seq=left_seq,
            right_aa_seq=right_seq,
        )
        Fasta2MFMIndex.create_and_save_index_from_fasta(
            fasta=hybridized_fasta, out_path=hybridized_fasta_fm_index
        )

        # Native run
        logger.info("Starting native Comet run")
        native_run = CometRun(
            mzml=mzml.path,
            fasta=hybrid_run_params.fasta,
            crux_comet_params=hybrid_run_params.crux_comet_params,
            out_dir=tmp_dir,
            scan_min=scan,
            scan_max=scan,
            num_threads=1,
        )
        native_run.run_comet_and_keep_only_results(crux_path=crux_path)
        move_file(
            src=native_run.standardized_comet_outputs.target,
            dest=native_out_dir / native_run.standardized_comet_outputs.target.name,
        )

        # # De-nativized Comet run
        # logger.info("Starting de-nativized Comet run")
        # denativized_run = CometRun(
        #     mzml=mzml.path,
        #     fasta=hybridized_fasta,
        #     crux_comet_params=hybrid_run_params.crux_comet_params,
        #     out_dir=tmp_dir,
        #     scan_min=scan,
        #     scan_max=scan,
        #     num_threads=1,
        # )
        # denativized_run.run_comet_and_keep_only_results(crux_path=crux_path)
        # move_file(
        #     src=denativized_run.standardized_comet_outputs.target,
        #     dest=expected_outputs[DENATIVIZED],
        # )

        # De-nativized HypedSearch run
        logger.info("Starting de-nativized HypedSearch run")
        denativized_params = deepcopy(hybrid_run_params)
        denativized_params.kmer_db_path = hybridized_kmer_db
        denativized_params.fasta = hybridized_fasta
        denativized_params.fasta_fm_index = hybridized_fasta_fm_index
        _, hybrid_run, _, _ = hybrid_run_on_spectrum(
            spectrum=mzml.get_spectrum(scan=scan),
            params=denativized_params,
            fasta_dir=tmp_dir,
            crux_path=crux_path,
            out_dir=tmp_dir,
        )
        move_file(
            src=hybrid_run.standardized_comet_outputs.target,
            dest=hybrid_out_dir / hybrid_run.standardized_comet_outputs.target.name,
        )


@click.command(
    name="run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    default=4,
    show_default=True,
    required=False,
    help="Number of cores to use to run experiment in parallel",
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the hybrid simulation config JSON",
)
@click.option(
    "--crux_path",
    "-cp",
    type=PathType(),
    required=False,
    help="Path to crux executable. If not provided, crux will be run via the Singularity container.",
)
@log_params
def cli_run_hybrid_finding_simulation_study_in_parallel(
    n_cores: int, config: Path, crux_path: Optional[Path]
):
    exp = HybridSimulationExperiment.load(path=config)
    exp.run_experiment_in_parallel(n_cores=n_cores, crux_path=crux_path)
    exp.combine_spectrum_txt_results()


@click.command(
    name="create-experiment-config",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--cut_method",
    "-cm",
    type=Union[str, float, int],
    required=True,
    help="",
)
@click.option(
    "--hybrid_run_params",
    "-hrp",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--psms",
    "-p",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--q_threshold",
    "-q",
    type=float,
    required=False,
    help="",
)
@log_params
def cli_create_experiment_config(
    hybrid_run_params: Path,
    psms: Path,
    cut_method: Union[str, float, int],
    q_threshold: Optional[float],
):
    pass


# @click.command(
#     name="run-in-serial",
#     context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
#     help="""
#     """,
# )
# @click.option(
#     "--config",
#     "-c",
#     type=PathType(),
#     required=True,
#     help="Path to the hybrid simulation config JSON",
# )
# @click.option(
#     "--crux_path",
#     "-cp",
#     type=PathType(),
#     required=False,
#     help="Path to crux executable. If not provided, crux will be run via the Singularity container.",
# )
# @log_params
# def cli_run_hybrid_finding_simulation_study_in_serial(
#     config: Path, crux_path: Optional[Path]
# ):
#     exp = HybridSimulationExperiment.load(path=config)
#     exp.run_experiment_in_serial(crux_path=crux_path)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    # cli.add_command(cli_run_hybrid_finding_simulation_study_in_serial)
    cli.add_command(cli_run_hybrid_finding_simulation_study_in_parallel)
    cli()
