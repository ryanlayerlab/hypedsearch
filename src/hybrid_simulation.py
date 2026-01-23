import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click

repo_dir = Path("/Users/erjo3868/repos/hypedsearch/hypedsearch")
os.chdir(repo_dir)
sys.path.append(str(repo_dir))
import random
import tempfile

from src.hypedsearch import (
    HybridFormer,
    HybridPSMScorer,
    HypedsearchRunConfig,
    SpectrumPreprocessor,
    hybrid_run_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide
from src.psm import CometPSM
from src.utils import setup_logger

logger = logging.getLogger(__name__)

DIR = Path(__file__).parent
random.seed(42069)


@dataclass
class NativeToHybridChange:
    """
    Class for representing turning a native sequence into a hybrid sequence within a protein.
    E.g., psm="AB" in prot="WXABYZ" -> hybrid="A-B" in "WAXYBZ"
    """

    left_hy_seq: str
    right_hy_seq: str
    prot_name: str
    original_prot_seq: str
    new_prot_seq: str

    @property
    def seq(self):
        return self.left_hy_seq + self.right_hy_seq


def make_native_seq_a_hybrid_seq_in_prot(
    native_seq: str, prot_seq: str
) -> Tuple[str, str, str]:
    # Split the hybrid sequence into two halves
    left_hy_seq, right_hy_seq = (
        native_seq[: int(len(native_seq) / 2)],
        native_seq[int(len(native_seq) / 2) :],
    )

    # Split the protein sequence into quarters and insert the hybrid halves between
    # the protein sequences first and second quarters and third and fourth quarters
    split_seq = prot_seq.split(native_seq)
    assert len(split_seq) == 2
    new_prot_seq = (
        split_seq[0][: int(len(split_seq[0]) / 2)]
        + left_hy_seq
        + split_seq[0][int(len(split_seq[0]) / 2) :]
        + split_seq[1][: int(len(split_seq[1]) / 2)]
        + right_hy_seq
        + split_seq[1][int(len(split_seq[1]) / 2) :]
    )
    return left_hy_seq, right_hy_seq, new_prot_seq


def validate_psm_for_hybrid_finding_simulation_study(
    psm: CometPSM,
    fasta: Fasta,
) -> bool:
    """
    Make sure the given PSM is easy to turn from a native to a hybrid which means:
    1) The PSM sequence appears in only one protein
    2) The PSM sequence appears only once in that protein
    If the PSM satisifies these conditions, return True, else False
    """

    # Make sure the PSM sequence appears in only one protein
    seq_containing_prots = list(fasta.proteins_that_contain_seqs([psm.seq])[psm.seq])
    if len(seq_containing_prots) != 1:
        return False

    # Make sure PSM appears only once in that protein
    prot_name = seq_containing_prots[0]
    if fasta.protein_name_to_seq_map[prot_name].count(psm.seq) != 1:
        return False

    return True


def prepare_psm_for_hybrid_finding_simulation_study(
    psm: CometPSM, fasta: Path
) -> Optional[NativeToHybridChange]:
    """
    Make sure the given PSM is easy to turn from a native to a hybrid which means:
    1) The PSM sequence appears in only one protein
    2) The PSM sequence appears only once in that protein
    If the PSM satisifies these conditions, turn it into a hybrid PSM
    """
    fasta = Fasta(path=fasta)
    if validate_psm_for_hybrid_finding_simulation_study(psm=psm, fasta=fasta):
        # Create the hybrid sequence by splitting the protein sequence at the PSM sequence location
        prot_name = list(fasta.proteins_that_contain_seqs([psm.seq])[psm.seq])[0]
        left_hy_seq, right_hy_seq, new_prot_seq = make_native_seq_a_hybrid_seq_in_prot(
            native_seq=psm.seq, prot_seq=fasta.protein_name_to_seq_map[prot_name]
        )
        return NativeToHybridChange(
            left_hy_seq=left_hy_seq,
            right_hy_seq=right_hy_seq,
            prot_name=prot_name,
            original_prot_seq=fasta.protein_name_to_seq_map[prot_name],
            new_prot_seq=new_prot_seq,
        )
    else:
        return None


def create_kmer_database_with_modified_protein(
    existing_kmer_db: Path, fasta: Path, modified_protein: Peptide, new_db_path: Path
) -> Tuple[List[Peptide], KmerDatabase]:
    kmer_db = KmerDatabase(db_path=existing_kmer_db)
    assert modified_protein.name in kmer_db.proteins
    kmer_db_prots = Fasta(path=fasta).get_proteins_by_name(
        names=list(
            set(prot for prot in kmer_db.proteins if prot != modified_protein.name)
        )
    ) + [modified_protein]
    new_kmer_db = KmerDatabase.create_db(
        db_path=new_db_path, proteins=kmer_db_prots, overwrite=True
    )
    return kmer_db_prots, new_kmer_db


def create_new_kmer_db_and_fasta_for_hybrid_simulation(
    native_to_hybrid_change: NativeToHybridChange,
    fasta: Path,
    existing_kmer_db: Path,
    new_kmer_db: Path,
    new_fasta: Path,
) -> KmerDatabase:
    kmer_db_prots, kmer_db = create_kmer_database_with_modified_protein(
        existing_kmer_db=existing_kmer_db,
        fasta=fasta,
        modified_protein=Peptide(
            seq=native_to_hybrid_change.new_prot_seq,
            name=native_to_hybrid_change.prot_name,
        ),
        new_db_path=new_kmer_db,
    )
    # Form hybrids
    Fasta.write_fasta(peptides=kmer_db_prots, path=new_fasta)
    return kmer_db


def create_hs_config_for_hybrid_simulation(
    hs_config_data: Dict,
    new_kmer_db: Path,
    out_dir: Optional[Path] = None,
) -> HypedsearchRunConfig:
    # Remove 'fasta' from 'psm_scorer' so there's no native-hybrid competition in hybrid run
    hs_config_data["psm_scorer"].pop("fasta", None)

    # Update kmer_db
    hs_config_data["hybrid_former"]["kmer_db"] = str(new_kmer_db)

    if out_dir is not None:
        hs_config_data["parent_out_dir"] = str(out_dir)

    return HypedsearchRunConfig(**hs_config_data)


def run_hybrid_simulation_on_psm(
    psm: CometPSM,
    spectrum: Spectrum,
    fasta: Path,
    out_dir: Path,
    kmer_db: Path,
    comet_params: Path,
):
    native_to_hybrid_change = prepare_psm_for_hybrid_finding_simulation_study(
        psm=psm, fasta=fasta
    )
    if native_to_hybrid_change is None:
        return None
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        new_kmer_db = tmp_dir / "kmers.db"
        new_fasta = tmp_dir / "proteins.fasta"
        _ = create_new_kmer_db_and_fasta_for_hybrid_simulation(
            native_to_hybrid_change=native_to_hybrid_change,
            fasta=fasta,
            existing_kmer_db=kmer_db,
            new_kmer_db=new_kmer_db,
            new_fasta=new_fasta,
        )
        _, comet_outputs = hybrid_run_on_spectrum(
            spectrum=spectrum,
            out_dir=out_dir,
            spectrum_preprocessor=SpectrumPreprocessor(),
            hybrid_former=HybridFormer(kmer_db=new_kmer_db, fasta=new_fasta),
            psm_scorer=HybridPSMScorer(comet_params=comet_params),
        )
    return comet_outputs


def processing_fcn(args):
    psm, spectrum, hs_config = args
    run_hybrid_simulation_on_psm(
        psm=psm,
        spectrum=spectrum,
        fasta=hs_config.fasta_path,
        out_dir=hs_config.hybrid_run_scan_results_dir,
        kmer_db=hs_config.kmer_db_path,
        comet_params=hs_config.comet_params_path,
    )


@click.command(
    name="run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
def cli_run_hybrid_finding_simulation_study(
    # sample_size: int
):
    hs_config_path = "results/hybrid_simulation/inputs/hs.config.json"
    hs_config = HypedsearchRunConfig.from_json(hs_config_path)
    mzml = Mzml(path=list(hs_config.mzml_to_scans.keys())[0])
    q_thresh = 0.01
    psms = CometPSM.from_txt(
        txt="results/1_1_Acet_Aspn_Islet_B35spike/native_run/assign-confidence.target.txt"
    )
    psms = [psm for psm in psms if psm.q_value <= q_thresh]
    logger.warning(f"There are {len(psms)} PSMs with q-value <= {q_thresh}.")

    # Sample N PSMs from confident ones
    sample_size = 500
    logger.warning(f"sample size = {sample_size}")
    top_psms = sorted(psms, key=lambda x: x.xcorr, reverse=True)[:sample_size]
    args = [(psm, mzml.id_to_spectrum[psm.spectrum_uid], hs_config) for psm in top_psms]
    logger.warning("Starting hybrid finding simulation study...")
    n_procs = 8
    with ProcessPoolExecutor(max_workers=n_procs) as exe:
        # list(exe.map(processing_fcn, args))
        futures = {exe.submit(processing_fcn, arg): arg for arg in args}
        total = len(futures)
        for i, fut in enumerate(as_completed(futures), start=1):
            arg = futures[fut]
            logger.warning(f"Finished PSM {i} of {total} (scan={arg[0].scan})")


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger(log_level=logging.WARNING)
    # cli.add_command(cli_create_native_run_snakemake_config)
    cli.add_command(cli_run_hybrid_finding_simulation_study)
    cli()
