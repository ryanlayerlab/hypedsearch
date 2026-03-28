import logging
import os
import random
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property, partial
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import click
from pydantic import BaseModel

from src.constants import DEFAULT_Q_THRESHOLD, MAC_CRUX_EXECUTABLE
from src.crux import CometRun, Crux
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Fasta2MFMIndex, Peptide
from src.psm import CometPSM
from src.utils import PathType, load_json, log_params, setup_logger, to_json

logger = logging.getLogger(__name__)

DIR = Path(__file__).parent
HALF = "half"
RANDOM = "random"
LEFT_SEQ = "left_hybrid_seq"
RIGHT_SEQ = "right_hybrid_seq"
NATIVE = "native"
DENATIVIZED = "denativized"
HYBRID = "hybrid"
CUT_METHODS = Literal[HALF, RANDOM]


def cut_seq_into_hybrid(seq: str, min_side_len: int, method: CUT_METHODS = HALF):
    assert len(seq) >= 2 * min_side_len
    if method == HALF:
        left_hy_seq, right_hy_seq = (
            seq[: int(len(seq) / 2)],
            seq[int(len(seq) / 2) :],
        )
    elif method == RANDOM:
        internal_seq = seq[min_side_len:-min_side_len]
        cut_idx = min_side_len + random.randint(0, len(internal_seq))
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
    aa_seq: str,
    left_aa_seq: str,
    right_aa_seq: str,
    hybridized_kmer_db: Union[str, Path],
    hybridized_fasta: Union[str, Path],
):
    assert (
        left_aa_seq + right_aa_seq == aa_seq
    ), f"Left and right sequences should concatenate to the original sequence. Got {left_aa_seq} and {right_aa_seq} which concatenate to {left_aa_seq + right_aa_seq}, not {aa_seq}."

    # Get the proteins that are in the existing kmer database. And for each of those
    # proteins remove the the given AA sequence from the proteins if it's present.
    kmer_db_prots = fasta.get_proteins_by_name(names=existing_kmer_db.proteins)
    new_prots = []
    num_prots_containing_seq = 0
    for pep in kmer_db_prots:
        new_pep = pep.model_copy(deep=True)
        if aa_seq in pep.seq:
            new_pep.seq = pep.seq.replace(aa_seq, "")
            num_prots_containing_seq += 1
        new_prots.append(new_pep)
    assert (
        num_prots_containing_seq > 0
    ), "The given sequence appears in none of the proteins in the kmer database! It's expected to be in at least one."

    # Add in left- and right- sequences as new proteins in the kmer database
    new_prots.extend(
        [
            Peptide(seq=left_aa_seq, name=LEFT_SEQ),
            Peptide(seq=right_aa_seq, name=RIGHT_SEQ),
        ]
    )

    # Create the new kmer database with the updated proteins
    KmerDatabase.create_db(
        db_path=hybridized_kmer_db, proteins=new_prots, overwrite=True
    )

    # Add all the other non-kmer databse proteins back into the FASTA file. Why?
    # Because when we run Comet, we want to run it on a FASTA
    new_prot_names = [prot.name for prot in new_prots]
    for prot in fasta.proteins:
        if prot.name not in new_prot_names:
            new_prots.append(prot)
    Fasta.write_fasta(peptides=new_prots, path=hybridized_fasta)


class HybridSimulator(BaseModel):
    cut_method: Literal[HALF, RANDOM] = HALF
    min_side_len: int = 3

    # Class and static methods
    @staticmethod
    def validate_aa_seq_for_hybrid_simulation(
        aa_seq: str,
        fasta: Fasta,
    ) -> bool:
        """
        Make sure the given PSM is easy to turn from a native to a hybrid which means:
        1) The PSM sequence appears in only one protein
        2) The PSM sequence appears only once in that protein
        If the PSM satisifies these conditions, return True, else False
        """
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

        return True

    # Instance methods
    def run_hybrid_simulation_on_spectrum(
        self,
        scan: int,
        aa_seq: str,
        mzml: Mzml,
        hybrid_run_params: HybridRunParams,
        out_dir: Path,
        crux_path: Optional[str | Path] = None,
    ):
        fasta = Fasta(path=hybrid_run_params.fasta)

        if not self.validate_aa_seq_for_hybrid_simulation(aa_seq=aa_seq, fasta=fasta):
            return
        left_hy_seq, right_hy_seq = cut_seq_into_hybrid(
            seq=aa_seq, method=self.cut_method, min_side_len=self.min_side_len
        )
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
                aa_seq=aa_seq,
                left_aa_seq=left_hy_seq,
                right_aa_seq=right_hy_seq,
            )
            Fasta2MFMIndex.create_and_save_index_from_fasta(
                fasta=hybridized_fasta, out_path=hybridized_fasta_fm_index
            )

            # Native Comet run
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
            # Move target output to output directory and add `native.` prefix
            txt = native_run.standardized_comet_outputs.target
            shutil.move(txt, out_dir / f"{NATIVE}.{txt.name}")

            # De-nativized Comet run
            denativized_run = CometRun(
                mzml=mzml.path,
                fasta=hybridized_fasta,
                crux_comet_params=hybrid_run_params.crux_comet_params,
                out_dir=tmp_dir,
                scan_min=scan,
                scan_max=scan,
                num_threads=1,
            )
            denativized_run.run_comet_and_keep_only_results(crux_path=crux_path)
            txt = denativized_run.standardized_comet_outputs.target
            shutil.move(txt, out_dir / f"{DENATIVIZED}.{txt.name}")

            # De-nativized HypedSearch run
            denativized_params = deepcopy(hybrid_run_params)
            denativized_params.kmer_db_path = hybridized_kmer_db
            denativized_params.fasta = hybridized_fasta
            denativized_params.fasta_fm_index = hybridized_fasta_fm_index
            _, hybrid_comet_run, _, _ = hybrid_run_on_spectrum(
                spectrum=mzml.get_spectrum(scan=scan),
                params=denativized_params,
                fasta_dir=tmp_dir,
                crux_path=crux_path,
                out_dir=tmp_dir,
            )
            txt = hybrid_comet_run.standardized_comet_outputs.target
            shutil.move(txt, out_dir / f"{HYBRID}.{txt.name}")


def run_hybrid_simulation_on_spectrum(
    scan: int,
    aa_seq: str,
    cut_method: Literal[HALF, RANDOM],
    mzml: Union[str, Path],
    hybrid_run_params: str | Path,
    out_dir: Union[str, Path],
    min_side_len: int = 3,
    crux_path: Optional[str | Path] = None,
):
    sim = HybridSimulator(
        cut_method=cut_method,
        min_side_len=min_side_len,
    )
    sim.run_hybrid_simulation_on_spectrum(
        aa_seq=aa_seq,
        scan=scan,
        mzml=Mzml(path=mzml),
        hybrid_run_params=HybridRunParams.load(path=hybrid_run_params),
        crux_path=crux_path,
        out_dir=out_dir,
    )


class HybridSimulationExperiment(BaseModel):
    hybrid_run_params: HybridRunParams
    hybrid_simulator: HybridSimulator
    scan_to_aa_seq: Dict[int, str]
    mzml: Path
    out_dir: Path

    @classmethod
    def load(cls, path: str | Path) -> "HybridSimulationExperiment":
        data = load_json(path=path)
        if isinstance(data["hybrid_run_params"], (str, Path)):
            data["hybrid_run_params"] = HybridRunParams.load(
                path=data["hybrid_run_params"]
            )
        data["hybrid_simulator"] = HybridSimulator(**data)
        return cls(**data)

    def save(self, path: str | Path):
        to_json(data=self.model_dump(mode="json"), path=path)

    def run_experiment(self, n_cores: int, crux_path: Optional[str | Path] = None):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        partial_fcn = partial(
            self.hybrid_simulator.run_hybrid_simulation_on_spectrum,
            mzml=Mzml(path=self.mzml),
            hybrid_run_params=self.hybrid_run_params,
            out_dir=self.out_dir,
            crux_path=crux_path,
        )
        with ProcessPoolExecutor(max_workers=n_cores) as ex:
            future_to_scan = {
                ex.submit(partial_fcn, scan, aa_seq): scan
                for scan, aa_seq in self.scan_to_aa_seq.items()
            }
            for idx, future in enumerate(as_completed(future_to_scan)):
                try:
                    _ = future.result()
                    logger.info(
                        f"Finished scan {future_to_scan[future]} ({idx + 1} of {len(self.scan_to_aa_seq)})"
                    )
                except Exception as e:
                    logger.warning(
                        f"Task failed for scan {future_to_scan[future]}: {e}"
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
def cli_run_hybrid_finding_simulation_study(
    n_cores: int, config: Path, crux_path: Optional[Path]
):
    exp = HybridSimulationExperiment.load(path=config)
    exp.run_experiment(n_cores=n_cores, crux_path=crux_path)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_run_hybrid_finding_simulation_study)
    cli()
