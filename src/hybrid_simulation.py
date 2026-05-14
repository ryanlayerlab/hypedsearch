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
from typing import ClassVar, Counter, Dict, List, Literal, Optional, Tuple, Union

import click
import numpy as np
import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import GIT_REPO_DIR, NATIVE
from src.crux import CometRun, Crux, run_comet
from src.hypedsearch import HybridRunParams, hybrid_run_on_spectrum
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Fasta2MFMIndex, Peptide
from src.psm import CometPSM, PeptideSeqSpectrumComparer
from src.utils import (
    CmdLineResult,
    PathType,
    flatten_list_of_lists,
    load_json,
    log_params,
    move_file,
    run_in_parallel,
    setup_logger,
    to_json,
)

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
    logger.debug("Creating hybridized k-mer database and FASTA for simulation")
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
        logger.debug(
            f"{validation_failure_msg} Sequence appears in multiple proteins: {seq_containing_prots}"
        )
        return False

    # Make sure sequence appears only once in that protein
    prot_name = seq_containing_prots[0]
    if fasta.protein_name_to_seq_map[prot_name].count(aa_seq) != 1:
        logger.debug(
            f"{validation_failure_msg} Sequence appears in one protein but appears > 1 times in {prot_name}."
        )
        return False

    # Make sure sequence each side length is >= min_side_len
    if len(left_seq) < min_side_len or len(right_seq) < min_side_len:
        logger.debug(
            f"{validation_failure_msg} One of the sides of the hybrid sequence is shorter than the minimum side length of {min_side_len}. Left seq: {left_seq}, right seq: {right_seq}."
        )
        return False
    return True


class HybridSimulationOutput(BaseModel):
    native_txt: Path
    denativized_txt: Path
    denativeized_hypedsearch_txt: Path


def cut_and_validate_seq(
    psm: CometPSM,
    cut_method: Union[RANDOM, float, int],
    min_side_len: int,
    fasta: Fasta,
) -> Tuple[int, str, str] | None:
    if len(psm.seq) < 2 * min_side_len:
        return None
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
        return (psm.scan, left_seq, right_seq)
    else:
        return None


class HybridSimulationExperiment(BaseModel):
    hybrid_run_params: HybridRunParams
    scan_to_left_right_seq: Dict[int, Tuple[str, str]]
    mzml: Mzml
    parent_out_dir: Path

    scan: ClassVar[str] = "scan"
    n_seq: ClassVar[str] = "native_seq"
    n_xcorr: ClassVar[str] = "native_xcorr"
    h_seq: ClassVar[str] = "hybrid_seq"
    h_xcorr: ClassVar[str] = "hybrid_xcorr"
    l_seq: ClassVar[str] = "left_seq"
    r_seq: ClassVar[str] = "right_seq"

    @field_validator("parent_out_dir", mode="before")
    @classmethod
    def create_output_dir(cls, v):
        path = Path(v)
        path.mkdir(parents=True, exist_ok=True)
        return path

    @cached_property
    def scan_to_left_seq(self) -> Dict[int, str]:
        return {
            scan: left_right_seq[0]
            for scan, left_right_seq in self.scan_to_left_right_seq.items()
        }

    @cached_property
    def scan_to_right_seq(self) -> Dict[int, str]:
        return {
            scan: left_right_seq[1]
            for scan, left_right_seq in self.scan_to_left_right_seq.items()
        }

    @cached_property
    def spectra(self) -> List[Spectrum]:
        return self.mzml.ms2_spectra

    @cached_property
    def scan_to_spectrum(self) -> Dict[int, Spectrum]:
        return {sp.scan: sp for sp in self.spectra}

    @staticmethod
    def create_native_out_dir(parent_out_dir: Path) -> Path:
        d = parent_out_dir / "native_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @staticmethod
    def create_hybrid_out_dir(parent_out_dir: Path) -> Path:
        d = parent_out_dir / "hybrid_run"
        d.mkdir(parents=True, exist_ok=True)
        return d

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
    def num_spectra(self):
        return len(self.scan_to_left_right_seq)

    @property
    def native_comet_txt(self) -> Path:
        native_txts = list(self.native_out_dir.glob("*.txt"))
        assert (
            len(native_txts) == 1
        ), f"Expected exactly one native txt file in {self.native_out_dir}, but found {len(native_txts)}"
        return native_txts[0]

    @property
    def hybrid_combined_txt(self) -> Path:
        return self.parent_out_dir / "combined_hybrid_psms.txt"

    @cached_property
    def native_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.native_comet_txt)

    @cached_property
    def hybrid_psms(self) -> List[CometPSM]:
        return CometPSM.from_txt(txt=self.hybrid_combined_txt)

    @cached_property
    def scan_to_sorted_native_psms(self) -> Dict[str, List[CometPSM]]:
        scan_to_psms = defaultdict(list)
        for psm in self.native_psms:
            scan_to_psms[psm.scan].append(psm)
        scan_to_psms = {
            scan: sorted(psms, key=lambda psm: psm.num)
            for scan, psms in scan_to_psms.items()
        }
        return scan_to_psms

    @cached_property
    def scan_to_sorted_hybrid_psms(self) -> Dict[str, List[CometPSM]]:
        scan_to_psms = defaultdict(list)
        for psm in self.hybrid_psms:
            scan_to_psms[psm.scan].append(psm)
        scan_to_psms = {
            scan: sorted(psms, key=lambda psm: psm.num)
            for scan, psms in scan_to_psms.items()
        }
        return scan_to_psms

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
        hybrid_run_params: str | Path | HybridRunParams,
        mzml: str | Path,
        psms: List[CometPSM],
        parent_out_dir: str | Path,
        cut_method: Union[RANDOM, float, int] = 0.5,
        min_side_len: int = 3,
        n_cores: int = 1,
    ):
        logger.info("Creating hybrid simulation experiment config")
        assert set(Counter(psm.scan for psm in psms).values()) == set(
            [1]
        ), "Each scan should have exactly one PSM associated with it"
        if isinstance(hybrid_run_params, (str, Path)):
            hybrid_run_params = HybridRunParams.load(path=hybrid_run_params)
        fasta = Fasta(path=hybrid_run_params.fasta)
        scan_to_left_right_seqs = {}
        validation_fcn = lambda psm: cut_and_validate_seq(
            psm=psm, cut_method=cut_method, min_side_len=min_side_len, fasta=fasta
        )
        results = run_in_parallel(
            fcn_of_one_variable=validation_fcn,
            input_array=psms,
            parallel_type="thread",
            n_cores=n_cores,
        )
        for result in results:
            if result is not None:
                scan, left_seq, right_seq = result
                scan_to_left_right_seqs[scan] = (left_seq, right_seq)

        # for idx, psm in enumerate(psms):
        #     print(f"Processing PSM {idx + 1} of {len(psms)}", end="\r")
        #     if len(psm.seq) < 2 * min_side_len:
        #         continue
        #     left_seq, right_seq = cut_seq_into_hybrid(
        #         seq=psm.seq,
        #         method=cut_method,
        #         min_side_len=min_side_len,
        #     )
        #     if validate_hybrid_for_hybrid_simulation(
        #         left_seq=left_seq,
        #         right_seq=right_seq,
        #         fasta=fasta,
        #         min_side_len=min_side_len,
        #     ):
        #         scan_to_left_right_seqs[psm.scan] = (left_seq, right_seq)
        #     else:
        #         logger.debug(
        #             f"PSM (scan={psm.scan}) with sequence {psm.seq} failed validation. Ignoring this PSM"
        #         )
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

    @property
    def shared_scans(self):
        return set(self.scan_to_sorted_native_psms.keys()).intersection(
            set(self.scan_to_sorted_hybrid_psms.keys())
        )

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
        data["mzml"] = str(self.mzml.path.relative_to(GIT_REPO_DIR))
        data["parent_out_dir"] = str(self.parent_out_dir.relative_to(GIT_REPO_DIR))
        to_json(
            data=data,
            path=path,
        )

    def run_hybrid_simulation_on_spectrum(
        self,
        scan: int,
        crux_path: Optional[str | Path] = None,
    ):
        run_hybrid_simulation_on_spectrum(
            scan=scan,
            left_seq=self.scan_to_left_seq[scan],
            right_seq=self.scan_to_right_seq[scan],
            mzml=self.mzml,
            hybrid_run_params=self.hybrid_run_params,
            out_dir=self.hybrid_out_dir,
            crux_path=crux_path,
        )

    def run_experiment_in_parallel(
        self, n_cores: int, crux_path: Optional[str | Path] = None
    ):
        partial_fcn = partial(
            self.run_hybrid_simulation_on_spectrum,
            crux_path=crux_path,
        )
        with ProcessPoolExecutor(max_workers=n_cores) as ex:
            future_to_scan = {
                ex.submit(partial_fcn, scan): scan
                for scan in self.scan_to_left_right_seq.keys()
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

    def combine_hybrid_txts(self):
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
        assert (
            self.scan_to_left_right_seq[scan][0] + self.scan_to_left_right_seq[scan][1]
            == native_psm.seq
        )
        return {
            self.scan: scan,
            self.n_seq: native_psm.seq,
            self.h_seq: hybrid_psm.seq,
            self.n_xcorr: native_psm.xcorr,
            self.h_xcorr: hybrid_psm.xcorr,
            self.l_seq: self.scan_to_left_right_seq[scan][0],
            self.r_seq: self.scan_to_left_right_seq[scan][1],
        }

    def summarize_results(self) -> Dict:
        results = defaultdict(set)
        for scan in self.shared_scans:
            # Get top PSMs
            native_psm = self.scan_to_sorted_native_psms[scan][0]
            assert (
                self.scan_to_left_seq[scan] + self.scan_to_right_seq[scan]
                == native_psm.seq
            )
            # Get the hybrid PSM with the same sequence if it exists
            matching_hybrid_psms = list(
                filter(
                    lambda psm: psm.seq == native_psm.seq,
                    self.scan_to_sorted_hybrid_psms[scan],
                )
            )
            if len(matching_hybrid_psms) == 0:
                results["missing"].add(scan)
            elif len(matching_hybrid_psms) == 1:
                hybrid_psm = matching_hybrid_psms[0]
                results[hybrid_psm.num].add(scan)
            else:
                raise ValueError(
                    f"Expected at most one hybrid PSM matching the native PSM sequence for scan {scan}, but found {len(matching_hybrid_psms)}"
                )
        return dict(results)

    @cached_property
    def top_native_vs_top_hybrid_df(self) -> pd.DataFrame:
        df = []
        for scan in self.shared_scans:
            # Get top PSMs
            native_psm = self.scan_to_sorted_native_psms[scan][0]
            hybrid_psm = self.scan_to_sorted_hybrid_psms[scan][0]
            assert (native_psm.num == 1) and (hybrid_psm.num == 1)
            df.append(
                self.compare_native_psm_to_hybrid_psm(
                    native_psm=native_psm,
                    hybrid_psm=hybrid_psm,
                )
            )
        df = pd.DataFrame(df)
        return df

    @property
    def native_equal_hybrid_mask(self) -> pd.Series:
        return (
            self.top_native_vs_top_hybrid_df[self.n_seq]
            == self.top_native_vs_top_hybrid_df[self.h_seq]
        )

    @property
    def num_spectra_where_top_hybrid_equals_top_native(self):
        return sum(self.native_equal_hybrid_mask)

    @property
    def native_not_equal_hybrid_df(self) -> pd.DataFrame:
        return self.top_native_vs_top_hybrid_df[~self.native_equal_hybrid_mask]

    @property
    def num_spectra_where_top_hybrid_beats_top_native(self):
        return self.native_not_equal_hybrid_df[
            self.native_not_equal_hybrid_df[self.h_xcorr]
            > self.native_not_equal_hybrid_df[self.n_xcorr]
        ].shape[0]

    def add_left_and_right_support_columns(
        self, df: pd.DataFrame, ppm_tol: float
    ) -> pd.DataFrame:
        expected_colms = [
            self.scan,
            self.n_seq,
            self.h_seq,
            self.n_xcorr,
            self.h_xcorr,
            self.l_seq,
            self.r_seq,
        ]
        assert set(expected_colms).issubset(
            set(df.columns)
        ), f"DataFrame should contain the following columns: {expected_colms}"
        scan_to_spectrum_native_psm_comparison = {
            row[self.scan]: PeptideSeqSpectrumComparer(
                spectrum=self.scan_to_spectrum[row[self.scan]],
                seq=row[self.n_seq],
                peak_to_ion_ppm_tol=ppm_tol,
            )
            for _, row in df.iterrows()
        }
        df["left_ion_support"] = df.apply(
            lambda row: scan_to_spectrum_native_psm_comparison[
                row[self.scan]
            ].ion_support_for_left_seq(left_seq=row[self.l_seq]),
            axis=1,
        )
        df["right_ion_support"] = df.apply(
            lambda row: scan_to_spectrum_native_psm_comparison[
                row[self.scan]
            ].ion_support_for_right_seq(right_seq=row[self.r_seq]),
            axis=1,
        )
        df["left_support"] = df.left_ion_support.apply(lambda ions: len(ions))
        df["right_support"] = df.right_ion_support.apply(lambda ions: len(ions))
        return df


def run_hybrid_simulation_experiment(
    hybrid_run_params: str | Path | HybridRunParams,
    mzml: str | Path,
    parent_out_dir: str | Path,
    n_cores: int,
    cut_method: Union[RANDOM, float, int] = 0.5,
    min_side_len: int = 3,
    crux_path: str | Path | None = None,
) -> HybridSimulationExperiment:
    if isinstance(hybrid_run_params, (str, Path)):
        hybrid_run_params = HybridRunParams.load(path=hybrid_run_params)
    # Run Comet natively
    logger.info("Starting native Comet run")
    native_comet_run = run_comet(
        fasta=hybrid_run_params.fasta,
        mzml=mzml,
        crux_comet_params=hybrid_run_params.crux_comet_params,
        out_dir=HybridSimulationExperiment.create_native_out_dir(
            parent_out_dir=parent_out_dir
        ),
        decoy_search=0,
        scan_min=0,
        scan_max=0,
        crux_path=crux_path,
    )
    psms = CometPSM.get_top_psms(
        psms=CometPSM.from_txt(txt=native_comet_run.standardized_comet_outputs.target)
    )

    # Create experiment config
    exp = HybridSimulationExperiment.create_experiment(
        hybrid_run_params=hybrid_run_params,
        mzml=mzml,
        psms=psms,
        parent_out_dir=parent_out_dir,
        cut_method=cut_method,
        min_side_len=min_side_len,
    )

    # Run experiment
    exp.run_experiment_in_parallel(n_cores=n_cores, crux_path=crux_path)

    # Combine hybrid txts
    exp.combine_hybrid_txts()
    return exp


def run_hybrid_simulation_on_spectrum(
    scan: int,
    left_seq: str,
    right_seq: str,
    mzml: Mzml,
    hybrid_run_params: HybridRunParams,
    out_dir: Path,
    crux_path: Optional[str | Path] = None,
) -> Path:
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

        # De-nativized HypedSearch run
        logger.debug(f"Starting de-nativized HypedSearch run for scan {scan}")
        denativized_params = deepcopy(hybrid_run_params)
        denativized_params.kmer_db_path = hybridized_kmer_db
        denativized_params.fasta = hybridized_fasta
        denativized_params.fasta_fm_index = hybridized_fasta_fm_index
        process, hybrid_comet_run, hybrid_seq_to_position_strs = hybrid_run_on_spectrum(
            spectrum=mzml.get_spectrum(scan=scan),
            params=denativized_params,
            fasta_dir=tmp_dir,
            crux_path=crux_path,
            out_dir=tmp_dir,
        )
        out_path = out_dir / hybrid_comet_run.standardized_comet_outputs.target.name
        move_file(
            src=hybrid_comet_run.standardized_comet_outputs.target,
            dest=out_path,
        )
    return out_path


@click.command(
    name="run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--num_cores",
    "-nc",
    type=int,
    required=True,
    help="",
)
@click.option(
    "--min_side_len",
    "-msl",
    type=int,
    required=False,
    default=3,
    show_default=True,
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
    "--crux_path",
    "-cp",
    type=PathType(),
    required=False,
    help="Path to crux executable. If not provided, crux will be run via the Singularity container.",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--cut_method",
    "-cm",
    type=str,
    required=True,
    default="half",
    show_default=True,
)
@log_params
def cli_run_hybrid_simulation(
    num_cores: int,
    hybrid_run_params: Path,
    crux_path: Path | None,
    mzml: Path,
    out_dir: Path,
    cut_method: str,
    min_side_len: int,
):
    if cut_method == "half":
        cut_method = 0.5
    elif cut_method == "random":
        cut_method = RANDOM
    else:
        raise ValueError(
            f"Invalid cut method: {cut_method}. Allowed values are 'half' and 'random'."
        )
    run_hybrid_simulation_experiment(
        hybrid_run_params=hybrid_run_params,
        mzml=mzml,
        parent_out_dir=out_dir,
        n_cores=num_cores,
        cut_method=cut_method,
        crux_path=crux_path,
        min_side_len=min_side_len,
    )


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
    cli.add_command(cli_run_hybrid_simulation)
    cli()
