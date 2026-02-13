import logging
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Literal, Tuple, Union

import click

from src.constants import DEFAULT_Q_THRESHOLD
from src.crux import Crux

repo_dir = Path("/Users/erjo3868/repos/hypedsearch/hypedsearch")
os.chdir(repo_dir)
sys.path.append(str(repo_dir))
import random
import tempfile

from src.hypedsearch import (
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide
from src.psm import CometPSM
from src.utils import setup_logger

logger = logging.getLogger(__name__)

DIR = Path(__file__).parent
HALF = "HALF"
RANDOM = "RANDOM"
CUT_METHODS = Literal[HALF, RANDOM]
random.seed(42069)


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


# @dataclass
# class NativeToHybridChange:
#     """
#     Class for representing turning a native sequence into a hybrid sequence within a protein.
#     E.g., psm="AB" in prot="WXABYZ" -> hybrid="A-B" in "WAXYBZ"
#     """

#     left_hy_seq: str
#     right_hy_seq: str
#     prot_name: str
#     original_prot_seq: str
#     new_prot_seq: str

#     @property
#     def seq(self):
#         return self.left_hy_seq + self.right_hy_seq


# def make_native_seq_a_hybrid_seq_in_prot(
#     native_seq: str, prot_seq: str
# ) -> Tuple[str, str, str]:
#     # Split the hybrid sequence into two halves
#     left_hy_seq, right_hy_seq = (
#         native_seq[: int(len(native_seq) / 2)],
#         native_seq[int(len(native_seq) / 2) :],
#     )

#     # Split the protein sequence into quarters and insert the hybrid halves between
#     # the protein sequences first and second quarters and third and fourth quarters
#     split_seq = prot_seq.split(native_seq)
#     assert len(split_seq) == 2
#     new_prot_seq = (
#         split_seq[0][: int(len(split_seq[0]) / 2)]
#         + left_hy_seq
#         + split_seq[0][int(len(split_seq[0]) / 2) :]
#         + split_seq[1][: int(len(split_seq[1]) / 2)]
#         + right_hy_seq
#         + split_seq[1][int(len(split_seq[1]) / 2) :]
#     )
#     return left_hy_seq, right_hy_seq, new_prot_seq


# def prepare_psm_for_hybrid_finding_simulation_study(
#     psm: CometPSM, fasta: Path
# ) -> Optional[NativeToHybridChange]:
#     """
#     Make sure the given PSM is easy to turn from a native to a hybrid which means:
#     1) The PSM sequence appears in only one protein
#     2) The PSM sequence appears only once in that protein
#     If the PSM satisifies these conditions, turn it into a hybrid PSM
#     """
#     fasta = Fasta(path=fasta)
#     if validate_psm_for_hybrid_finding_simulation_study(psm=psm, fasta=fasta):
#         # Create the hybrid sequence by splitting the protein sequence at the PSM sequence location
#         prot_name = list(fasta.proteins_that_contain_seqs([psm.seq])[psm.seq])[0]
#         left_hy_seq, right_hy_seq, new_prot_seq = make_native_seq_a_hybrid_seq_in_prot(
#             native_seq=psm.seq, prot_seq=fasta.protein_name_to_seq_map[prot_name]
#         )
#         return NativeToHybridChange(
#             left_hy_seq=left_hy_seq,
#             right_hy_seq=right_hy_seq,
#             prot_name=prot_name,
#             original_prot_seq=fasta.protein_name_to_seq_map[prot_name],
#             new_prot_seq=new_prot_seq,
#         )
#     else:
#         return None


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
    validation_failure_msg = f"Validation failed for PSM {psm.uid}"
    seq_containing_prots = list(fasta.proteins_that_contain_seqs([psm.seq])[psm.seq])
    if len(seq_containing_prots) != 1:
        logger.info(validation_failure_msg)
        return False

    # Make sure PSM appears only once in that protein
    prot_name = seq_containing_prots[0]
    if fasta.protein_name_to_seq_map[prot_name].count(psm.seq) != 1:
        logger.info(validation_failure_msg)
        return False

    return True


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


def create_new_kmer_db_and_fasta_for_simulation_with_hybrid_within_prot(
    existing_kmer_db: Union[str, Path],
    new_kmer_db: Union[str, Path],
    fasta: Union[str, Path, Fasta],
    psm_seq: str,
    left_seq: str,
    right_seq: str,
    new_fasta: Union[str, Path],
):
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
    existing_kmer_db: Union[str, Path],
    new_kmer_db: Union[str, Path],
    fasta: Union[str, Path, Fasta],
    psm_seq: str,
    left_seq: str,
    right_seq: str,
    new_fasta: Union[str, Path],
):
    if isinstance(fasta, (str, Path)):
        fasta = Fasta(path=fasta)
    kmer_db = KmerDatabase(db_path=existing_kmer_db)
    kmer_db_prots = fasta.get_proteins_by_name(names=kmer_db.proteins)
    new_prots = []
    num_prots_containing_seq = 0
    for pep in kmer_db_prots:
        new_pep = pep.model_copy(deep=True)
        if psm_seq in pep.seq:
            new_pep.seq = pep.seq.replace(psm_seq, "")
            num_prots_containing_seq += 1
        new_prots.append(new_pep)
    assert num_prots_containing_seq > 0
    logger.info(
        f"Number of proteins containing seq {psm_seq}: {num_prots_containing_seq}"
    )
    # Add in prot as hybrid
    new_prots.extend(
        [
            Peptide(seq=left_seq, name="left_hybrid_seq"),
            Peptide(seq=right_seq, name="right_hybrid_seq"),
        ]
    )
    KmerDatabase.create_db(db_path=new_kmer_db, proteins=new_prots, overwrite=True)

    new_prot_names = [prot.name for prot in new_prots]
    for prot in fasta.proteins:
        if prot.name not in new_prot_names:
            new_prots.append(prot)
    Fasta.write_fasta(peptides=new_prots, path=new_fasta)


@dataclass
class HybridSimulation:
    config: HypedsearchRunConfig
    assign_confidence_txt: Path
    parent_out_dir: Path
    q_threshold: int = DEFAULT_Q_THRESHOLD

    def __post_init__(self):
        if isinstance(self.config, (str, Path)):
            self.config = HypedsearchRunConfig.from_json(self.config)
        self.parent_out_dir = Path(self.parent_out_dir)
        self.parent_out_dir.mkdir(parents=True, exist_ok=True)
        if self.missing_top_peptide_dir.exists():
            shutil.rmtree(self.missing_top_peptide_dir)
        self.missing_top_peptide_dir.mkdir(parents=True, exist_ok=True)
        if self.native_dir.exists():
            shutil.rmtree(self.native_dir)
        self.native_dir.mkdir(parents=True, exist_ok=True)
        if self.hybrid_dir.exists():
            shutil.rmtree(self.hybrid_dir)
        self.hybrid_dir.mkdir(parents=True, exist_ok=True)

    @property
    def missing_top_peptide_dir(self):
        return self.parent_out_dir / "missing_top_peptide"

    @property
    def native_dir(self):
        return self.parent_out_dir / "native_run"

    @property
    def hybrid_dir(self):
        return self.parent_out_dir / "hybrid_run"

    @cached_property
    def confident_native_psms(self):
        return sorted(
            [
                psm
                for psm in CometPSM.from_txt(txt=self.assign_confidence_txt)
                if psm.q_value <= self.q_threshold
            ],
            key=lambda psm: psm.xcorr,
            reverse=True,
        )

    def run_hybrid_simulation_on_psm(
        self,
        psm: CometPSM,
        spectrum: Spectrum,
        mzml_path: Union[str, Path],
        cut_method: CUT_METHODS = "HALF",
        min_side_len: int = 3,
    ):
        mzml = Mzml(path=mzml_path)
        fasta = Fasta(path=self.config.fasta)
        assert psm.sample == mzml.name, "PSM sample and MZML name should be the same!"
        if validate_psm_for_hybrid_finding_simulation_study(psm=psm, fasta=fasta):
            left_hy_seq, right_hy_seq = cut_seq_into_hybrid(
                seq=psm.seq, method=cut_method, min_side_len=min_side_len
            )
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_dir = Path(tmp_dir)
                new_kmer_db = tmp_dir / "kmers.db"
                new_fasta = tmp_dir / "proteins.fasta"
                # create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots(
                #     existing_kmer_db=self.config.kmer_db,
                #     new_kmer_db=new_kmer_db,
                #     new_fasta=new_fasta,
                #     fasta=self.config.fasta,
                #     psm_seq=psm.seq,
                #     left_seq=left_hy_seq,
                #     right_seq=right_hy_seq,
                # )
                create_new_kmer_db_and_fasta_for_simulation_with_hybrid_within_prot(
                    existing_kmer_db=self.config.kmer_db,
                    new_kmer_db=new_kmer_db,
                    new_fasta=new_fasta,
                    fasta=self.config.fasta,
                    psm_seq=psm.seq,
                    left_seq=left_hy_seq,
                    right_seq=right_hy_seq,
                )
                # Run HypedSearch
                new_config = self.config.model_copy(deep=True)
                new_config.kmer_db = new_kmer_db
                new_config.fasta = new_fasta
                hybrid_run_on_spectrum(
                    spectrum=spectrum,
                    hybrid_former=new_config.hybrid_former,
                    psm_scorer=new_config.psm_scorer,
                    out_dir=self.hybrid_dir,
                )

                # Run Comet on native FASTA and hybridized FASTA
                Crux().run_comet(
                    mzml=mzml_path,
                    fasta=self.config.fasta,
                    crux_comet_params=self.config.crux_comet_params,
                    out_dir=self.native_dir,
                    file_root=mzml.name,
                    scan_min=psm.scan,
                    scan_max=psm.scan,
                    num_threads=1,
                )
                Crux().run_comet(
                    mzml=mzml_path,
                    fasta=new_fasta,
                    crux_comet_params=self.config.crux_comet_params,
                    out_dir=self.missing_top_peptide_dir,
                    file_root=mzml.name,
                    scan_min=psm.scan,
                    scan_max=psm.scan,
                    num_threads=1,
                )


def processing_fcn(
    sim: HybridSimulation, psm: CometPSM, spectrum: Spectrum, mzml_path: Path
):
    # sim = HybridSimulation(
    #     config="results/tutorial/inputs/hs.config.json",
    #     assign_confidence_txt="results/tutorial/native_run/tutorial-assign-confidence.txt",
    #     parent_out_dir="results/hybrid_simulation",
    # )
    sim.run_hybrid_simulation_on_psm(
        psm=psm,
        spectrum=spectrum,
        mzml_path=mzml_path,
    )


@click.command(
    name="run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--sample_size",
    "-n",
    type=int,
    default=100,
    show_default=True,
    required=False,
    help="",
)
def cli_run_hybrid_finding_simulation_study(sample_size: int):
    mzml_path = Path("data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML")
    uid_to_spectrum = Mzml(path=mzml_path).id_to_spectrum
    sim = HybridSimulation(
        config="results/tutorial/inputs/hs.config.json",
        assign_confidence_txt="results/tutorial/native_run/tutorial-assign-confidence.txt",
        parent_out_dir="results/020226_hybrid_simulation_hybridize_within_protein",
    )
    args = [
        (sim, psm, uid_to_spectrum[psm.spectrum_uid], mzml_path)
        for psm in sim.confident_native_psms[:sample_size]
    ]
    logger.warning("Starting hybrid finding simulation study...")
    n_procs = 8
    with ProcessPoolExecutor(max_workers=n_procs) as exe:
        # list(exe.map(processing_fcn, args))
        futures = {exe.submit(processing_fcn, *arg): arg for arg in args}
        total = len(futures)
        for i, fut in enumerate(as_completed(futures), start=1):
            arg = futures[fut]
            logger.warning(f"Finished PSM {i} of {total} (scan={arg[1].scan})")


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    # cli.add_command(cli_create_native_run_snakemake_config)
    cli.add_command(cli_run_hybrid_finding_simulation_study)
    cli()
