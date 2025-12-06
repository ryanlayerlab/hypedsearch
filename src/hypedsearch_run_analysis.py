import itertools
import time
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from functools import cached_property
from pathlib import Path
from typing import ClassVar, Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from pydantic import BaseModel, Field, model_validator
from scipy.interpolate import PchipInterpolator
from scipy.stats import ecdf
from statsmodels.distributions.empirical_distribution import ECDF

from src.comet_utils import CometPSM
from src.constants import (
    DECOY,
    DEFAULT_FPR,
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_RANGE,
    DEFAULT_SCORE_CHANGE_RANGE,
    HYBRID,
    NATIVE,
    Q_VAL,
    Q_VAL_THRESH,
    SPECTRA_PSMS_FILE_NAME,
    TARGET,
    XCORR,
)
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import HypedsearchRunConfig
from src.kmer_database import KmerDatabase
from src.mass_spectra import Spectrum
from src.peptide_spectrum_comparison import PSM
from src.peptides_and_ions import Fasta
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.protein_abundance import ProteinAbundance
from src.utils import (
    PathType,
    Position,
    decompress_and_depickle,
    flatten_list_of_lists,
    get_positions_of_subseq_in_seq,
    pickle_and_compress,
    setup_logger,
)

SCORES = [
    "xcorr",
    "prop_intensity_supported",
    "prop_ions_matched",
    "prop_prefixes_supported",
    "prop_suffixes_supported",
    "q_value",
]


class HybridPSM(BaseModel):
    hybrids: List[HybridPeptide]
    psm: PSM


class RunResultsDir(BaseModel):
    folder: Path
    target_txts: List[Path]
    decoy_txts: List[Path]
    assign_conf_txt: Optional[Path] = None

    @classmethod
    def from_path(cls, path: Path) -> "RunResultsDir":
        path = Path(path)
        target_txts = []
        decoy_txts = []
        assign_conf_txt = None
        for txt in path.glob("*.txt"):
            if "assign-confidence" in txt.name:
                assert (
                    assign_conf_txt is None
                ), f"Multiple assign-confidence txts found in {path}"
                assign_conf_txt = txt
            elif TARGET in txt.name:
                target_txts.append(txt)
            elif DECOY in txt.name:
                decoy_txts.append(txt)
            else:
                raise RuntimeError(
                    f"Found a .txt that I don't know what to do with: {txt.name}"
                )
        return cls(
            folder=path,
            target_txts=target_txts,
            decoy_txts=decoy_txts,
            assign_conf_txt=assign_conf_txt,
        )


class SpectrumCometResults(BaseModel):
    targets: List[CometPSM] = Field(default_factory=list)
    decoys: List[CometPSM] = Field(default_factory=list)
    assign_conf: Optional[CometPSM] = None

    @model_validator(mode="after")
    def check_consistent_spectrum_uid(self):
        # Collect UIDs from targets and decoys
        uids = {psm.spectrum_uid for psm in self.targets + self.decoys}
        if len(uids) != 1:
            raise ValueError(f"Inconsistent spectrum_uid values found: {uids}")
        return self

    @property
    def spectrum_uid(self) -> str:
        if len(self.targets) > 0:
            return self.targets[0].spectrum_uid
        elif len(self.decoys) > 0:
            return self.decoys[0].spectrum_uid
        else:
            raise ValueError("No targets or decoys to get spectrum_uid from.")

    @cached_property
    def target(self) -> Optional[CometPSM]:
        if len(self.targets) == 0:
            return None
        else:
            top_psms = CometPSM.get_top_psms(psms=self.targets)
            assert (
                len(top_psms) == 1
            ), f"Expected 1 top target PSM, found {len(top_psms)}"
            return top_psms[0]

    @cached_property
    def decoy(self) -> CometPSM:
        if len(self.decoys) == 0:
            return None
        else:
            top_psms = CometPSM.get_top_psms(psms=self.decoys)
            assert (
                len(top_psms) == 1
            ), f"Expected 1 top decoy PSM, found {len(top_psms)}"
            return top_psms[0]

    @classmethod
    def from_dir(
        cls, path: Path, by_spectrum: bool = True
    ) -> Union[List["SpectrumCometResults"], Dict[str, "SpectrumCometResults"]]:
        # Load targets
        target_txts = []
        decoy_txts = []
        assign_conf_txt = None
        for txt in path.glob("*.txt"):
            if "assign-confidence" in txt.name:
                assert (
                    assign_conf_txt is None
                ), f"Multiple assign-confidence txts found in {path}"
                assign_conf_txt = txt
            elif TARGET in txt.name:
                target_txts.append(txt)
            elif DECOY in txt.name:
                decoy_txts.append(txt)
            else:
                raise RuntimeError(
                    f"Found a .txt that I don't know what to do with: {txt.name}"
                )

        spectrum_to_targets = CometPSM.from_txts(txts=target_txts, by_spectrum=True)
        spectrum_to_decoys = CometPSM.from_txts(txts=decoy_txts, by_spectrum=True)
        spectrum_to_assign_conf = (
            Spectrum.organize_by_spectrum(data=CometPSM.from_txt(txt=assign_conf_txt))
            if assign_conf_txt is not None
            else {}
        )
        results = []
        for spectrum_uid in set(
            list(spectrum_to_targets.keys())
            + list(spectrum_to_decoys.keys())
            + list(spectrum_to_assign_conf.keys())
        ):
            targets = spectrum_to_targets.get(spectrum_uid, [])
            decoys = spectrum_to_decoys.get(spectrum_uid, [])
            assign_conf = spectrum_to_assign_conf.get(spectrum_uid, None)
            results.append(cls(targets=targets, decoys=decoys, assign_conf=assign_conf))

        if by_spectrum:
            return Spectrum.organize_by_spectrum(data=results)
        else:
            return results

    def target_psm(
        self,
        spectrum: Spectrum,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> PSM:
        if self.assign_conf is not None:
            return PSM.from_spectrum_and_comet_psm(
                spectrum=spectrum,
                comet_psm=self.assign_conf,
                peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            )
        else:
            return PSM.from_spectrum_and_comet_psm(
                spectrum=spectrum,
                comet_psm=self.target,
                peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            )

    def decoy_psm(
        self,
        spectrum: Spectrum,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> PSM:
        return PSM.from_spectrum_and_comet_psm(
            spectrum=spectrum,
            comet_psm=self.decoy,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
        )


@dataclass
class ProteinRange(Position):
    protein: str

    @classmethod
    def from_pos(cls, protein: str, pos: Position):
        return cls(
            protein=protein,
            inclusive_start=pos.inclusive_start,
            exclusive_end=pos.exclusive_end,
        )


@dataclass
class HybridPosition:
    left: ProteinRange
    right: ProteinRange


class SpectrumPSMs(BaseModel):
    spectrum: Union[str, Spectrum]
    native_target: Optional[PSM] = None
    native_decoy: Optional[PSM] = None
    hybrid_target: Optional[PSM] = None
    hybrid_decoy: Optional[PSM] = None
    hybrids: Optional[List[HybridPeptide]] = None

    @property
    def spectrum_uid(self) -> str:
        if isinstance(self.spectrum, str):
            return self.spectrum
        else:
            return self.spectrum.uid

    @property
    def native_seq(self):
        return self.native_target.seq

    @property
    def native_q(self):
        if self.native_target is None:
            return None
        return self.native_target.q_value

    @property
    def native_xcorr(self):
        return self.native_target.xcorr

    @property
    def hybrid_seq(self):
        if self.hybrid_target is not None:
            return self.hybrid_target.seq
        else:
            return None

    @property
    def hybrid_xcorr(self):
        if self.hybrid_target is not None:
            return self.hybrid_target.xcorr
        else:
            return None

    @property
    def hybrid_hyphen_seqs(self):
        if self.hybrids is not None:
            return set(hy.hyphen_seq for hy in self.hybrids)
        else:
            return None

    @property
    def mzml(self):
        return self.spectrum.mzml

    def to_dicts(self) -> List[Dict]:
        row = {
            "spectrum": self.spectrum_uid,
            "precursor_mz": self.spectrum.precursor_mz,
            "precursor_charge": self.spectrum.precursor_charge,
            "retention_time": self.spectrum.retention_time,
        }
        for psm_type in [
            "native_target",
            "native_decoy",
            "hybrid_target",
            "hybrid_decoy",
        ]:
            for attr in SCORES + ["seq"]:
                psm = getattr(self, psm_type)
                row[f"{psm_type}_{attr}"] = (
                    getattr(psm, attr) if psm is not None else None
                )

        row["hybrid_seq"] = None
        row["left_prots"] = None
        row["right_prots"] = None
        if self.hybrids:
            rows = []
            for hy in self.hybrids:
                new_row = deepcopy(row)
                new_row["hybrid_seq"] = hy.hyphen_seq
                new_row["left_prots"] = ";".join(hy.left_proteins)
                new_row["right_prots"] = ";".join(hy.right_proteins)
                rows.append(new_row)
            return rows

        else:
            return [row]

    @staticmethod
    def to_df(spectra_psms: List["SpectrumPSMs"]) -> pd.DataFrame:
        data = flatten_list_of_lists(
            [spectrum_psms.to_dicts() for spectrum_psms in spectra_psms]
        )
        return pd.DataFrame(data)

    @staticmethod
    def save(
        spectra_psms: List["SpectrumPSMs"],
        out_path: Union[str, Path],
        overwrite: bool = False,
    ):
        out_path = Path(out_path)
        if out_path.exists() and not overwrite:
            logger.info(
                f"Spectra PSMs file already exists at {out_path}, not overwriting. "
                + "Skipping saving..."
            )
        else:
            data = [spectrum_psms.model_dump() for spectrum_psms in spectra_psms]
            out_path.parent.mkdir(parents=True, exist_ok=True)
            pickle_and_compress(obj=data, path=out_path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> Dict[str, "SpectrumPSMs"]:
        data = decompress_and_depickle(path=path)
        data = [cls(**d) for d in data]
        # return {d.spectrum_uid: d for d in data}
        return data

    @classmethod
    def from_native_and_hybrid_dirs(
        cls,
        native_dir: Path,
        hybrid_dir: Path,
        spectra: Dict[str, Spectrum],
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> List["SpectrumPSMs"]:
        native_spectra_results = SpectrumCometResults.from_dir(
            path=native_dir, by_spectrum=True
        )
        hybrid_spectra_results = SpectrumCometResults.from_dir(
            path=hybrid_dir, by_spectrum=True
        )
        all_spectra_uids = set(
            list(native_spectra_results.keys()) + list(hybrid_spectra_results.keys())
        )
        data = []
        logger.info(
            f"Creating SpectrumPSM objects for each spectrum (n={len(all_spectra_uids)})"
        )
        for spectrum_uid in all_spectra_uids:
            spectrum = spectra[spectrum_uid]
            native_target, native_decoy = None, None
            if spectrum_uid in native_spectra_results:
                native_target = native_spectra_results[spectrum_uid].target_psm(
                    spectrum=spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
                native_decoy = native_spectra_results[spectrum_uid].decoy_psm(
                    spectrum=spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
            hybrid_target, hybrid_decoy, hybrids = None, None, None
            if spectrum_uid in hybrid_spectra_results:
                if hybrid_spectra_results[spectrum_uid].target.is_hybrid:
                    hybrid_target = hybrid_spectra_results[spectrum_uid].target_psm(
                        spectrum=spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                    )
                hybrid_decoy = hybrid_spectra_results[spectrum_uid].decoy_psm(
                    spectrum=spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
            data.append(
                cls(
                    spectrum=spectrum,
                    native_target=native_target,
                    native_decoy=native_decoy,
                    hybrid_target=hybrid_target,
                    hybrid_decoy=hybrid_decoy,
                )
            )
        return data

    @classmethod
    def from_comet_txts(
        cls,
        native_target_txts: List[Path],
        native_decoy_txts: List[Path],
        native_assign_confidence_txt: Path,
        hybrid_target_txts: List[Path],
        hybrid_decoy_txts: List[Path],
        hybrid_assign_confidence_txt: Path,
        spectra: List[Spectrum],
        kmer_to_proteins_map: Dict[str, List[str]],
        peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
        min_side_len: int = DEFAULT_MIN_SIDE_LEN,
    ):
        # Load Comet PSMs
        logger.info("Loading Comet PSMs from txt files...")
        native_target_psms = cls.load_num1_comet_psms_from_txts(txts=native_target_txts)
        native_decoy_psms = cls.load_num1_comet_psms_from_txts(txts=native_decoy_txts)
        native_assign_confidence_psms = cls.load_num1_comet_psms_from_txts(
            txts=[native_assign_confidence_txt]
        )
        hybrid_target_psms = cls.load_num1_comet_psms_from_txts(
            txts=hybrid_target_txts, only_include_hybrids=True
        )
        hybrid_decoy_psms = cls.load_num1_comet_psms_from_txts(txts=hybrid_decoy_txts)
        hybrid_assign_confidence_psms = cls.load_num1_comet_psms_from_txts(
            txts=[hybrid_assign_confidence_txt], only_include_hybrids=True
        )

        # Create SpectrumPSMs objects
        all_spectra_ids = set(
            list(native_target_psms.keys())
            + list(native_decoy_psms.keys())
            + list(native_assign_confidence_psms.keys())
            + list(hybrid_target_psms.keys())
            + list(hybrid_decoy_psms.keys())
            + list(hybrid_assign_confidence_psms.keys())
        )
        spectra = {
            spectrum.uid: spectrum
            for spectrum in spectra
            if spectrum.uid in all_spectra_ids
        }
        logger.info(
            f"Creating SpectrumPSM objects for each spectrum (n={len(spectra)})"
        )
        psms = []
        for idx, val in enumerate(spectra.items()):
            if idx % 500 == 0:
                logger.info(f"Processing spectrum {idx + 1}/{len(spectra)}")
            spectrum_id, spectrum = val
            # Get hybrids for hybrid target PSM
            hybrid_target = None
            hybrid_assign_confidence = None
            hybrids = None
            if spectrum_id in hybrid_target_psms:
                comet_psm = hybrid_target_psms[spectrum_id]
                all_hybrids = find_possible_hybrids(
                    seq=comet_psm.seq,
                    kmer_to_proteins_map=kmer_to_proteins_map,
                    min_side_len=min_side_len,
                )
                if len(all_hybrids) > 0:
                    hybrid_target = PSM.from_spectrum_and_comet_psm(
                        comet_psm=comet_psm,
                        spectrum=spectrum,
                        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                    )
                    hybrid_assign_confidence = (
                        PSM.from_spectrum_and_comet_psm(
                            spectrum=spectrum,
                            comet_psm=hybrid_assign_confidence_psms[spectrum_id],
                            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                        )
                        if spectrum_id in hybrid_assign_confidence_psms
                        else None
                    )
                    HybridPeptide.set_proteins(
                        hybrids=all_hybrids, seq_to_proteins=kmer_to_proteins_map
                    )
                    hybrids = all_hybrids

            psm = cls(
                spectrum=spectrum,
                hybrids=hybrids,
                native_target=(
                    PSM.from_spectrum_and_comet_psm(
                        spectrum=spectrum,
                        comet_psm=native_target_psms[spectrum_id],
                        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                    )
                    if spectrum_id in native_target_psms
                    else None
                ),
                native_decoy=(
                    PSM.from_spectrum_and_comet_psm(
                        spectrum=spectrum,
                        comet_psm=native_decoy_psms[spectrum_id],
                        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                    )
                    if spectrum_id in native_decoy_psms
                    else None
                ),
                native_assign_confidence=(
                    PSM.from_spectrum_and_comet_psm(
                        spectrum=spectrum,
                        comet_psm=native_assign_confidence_psms[spectrum_id],
                        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                    )
                    if spectrum_id in native_assign_confidence_psms
                    else None
                ),
                hybrid_target=hybrid_target,
                hybrid_decoy=(
                    PSM.from_spectrum_and_comet_psm(
                        spectrum=spectrum,
                        comet_psm=hybrid_decoy_psms[spectrum_id],
                        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
                    )
                    if spectrum_id in hybrid_decoy_psms
                    else None
                ),
                hybrid_assign_confidence=hybrid_assign_confidence,
            )
            psms.append(psm)
        return psms

    @staticmethod
    def add_hybrids_to_spectra_psms(
        spectra_psms: List["SpectrumPSMs"],
        db_path: Path,
        min_side_len: int = DEFAULT_MIN_SIDE_LEN,
    ) -> List["SpectrumPSMs"]:
        kmer_to_proteins_map = kmer_to_proteins_map = KmerDatabase(
            db_path=db_path
        ).kmer_to_proteins_map.kmer_to_protein_map
        for psms in spectra_psms:
            if psms.hybrid_target is not None:
                psms.hybrid_target.seq
                hybrids = find_possible_hybrids(
                    seq=psms.hybrid_target.seq,
                    kmer_to_proteins_map=kmer_to_proteins_map,
                    min_side_len=min_side_len,
                )
                if len(hybrids) > 0:
                    psms.hybrids = hybrids
        return spectra_psms

    @classmethod
    def from_hs_config(
        cls,
        hs_config: Union[HypedsearchRunConfig, Path],
        min_side_len: int = DEFAULT_MIN_SIDE_LEN,
        save: bool = False,
    ) -> List["SpectrumPSMs"]:
        if isinstance(hs_config, Path):
            hs_config = HypedsearchRunConfig.from_json(path=hs_config)

        logger.info("Loading spectra")
        spectra = Spectrum.organize_by_spectrum(
            Spectrum.load_spectra_from_mzmls(
                mzmls=list(hs_config.mzml_to_scans.keys()),
            )
        )
        # Form SpectrumPSMs objects from native and hybrid dirs
        start_time = time.perf_counter()
        spectra_psms = SpectrumPSMs.from_native_and_hybrid_dirs(
            native_dir=hs_config.native_run_dir,
            hybrid_dir=hs_config.hybrid_run_dir,
            spectra=spectra,
            peak_to_ion_ppm_tol=hs_config.hybrid_former.peak_to_ion_ppm_tol,
        )
        dur = time.perf_counter() - start_time
        logger.info(
            f"Creating SpectrumPSM objects for {len(spectra_psms)} spectra took {dur:.2f} seconds"
        )
        logger.info("Adding hybrids to SpectrumPSMs objects")
        start_time = time.perf_counter()
        spectra_psms = cls.add_hybrids_to_spectra_psms(
            spectra_psms=spectra_psms,
            db_path=hs_config.hybrid_former.kmer_db,
            min_side_len=min_side_len,
        )
        dur = time.perf_counter() - start_time
        logger.info(f"Adding hybrids took {dur:.2f} seconds")

        if save:
            out_path = SpectrumPSMs.default_save_path(
                out_dir=hs_config.parent_out_dir, min_side_len=min_side_len
            )
            logger.info(f"Saving spectra to {out_path}")
            start_time = time.perf_counter()
            SpectrumPSMs.save(
                spectra_psms=spectra_psms, out_path=out_path, overwrite=True
            )
            dur = time.perf_counter() - start_time
            logger.info(f"Saving {len(spectra_psms)} spectra took {dur:.2f} seconds")
        return spectra_psms

    @staticmethod
    def default_save_path(out_dir: Path, min_side_len: int):
        return out_dir / f"spectra_psms_minSideLen{min_side_len}.pklz"

    @staticmethod
    def filter_to_has_hybrid(
        psms: List["SpectrumPSMs"],
    ) -> List["SpectrumPSMs"]:
        return list(
            filter(
                lambda psm: (psm.hybrid_target is not None)
                and (psm.hybrids is not None),
                psms,
            )
        )

    @staticmethod
    def filter_to_has_native_q(
        psms: List["SpectrumPSMs"],
    ) -> List["SpectrumPSMs"]:
        return [psm for psm in psms if psm.native_q is not None]

    @staticmethod
    def get_winners_by_xcorr(
        psms: List["SpectrumPSMs"], psm_type: Literal[NATIVE, HYBRID]
    ) -> List["SpectrumPSMs"]:
        psms = SpectrumPSMs.filter_to_has_hybrid(psms=psms)
        if psm_type == NATIVE:
            return [psm for psm in psms if (psm.native_xcorr > psm.hybrid_xcorr)]
        elif psm_type == HYBRID:
            return [psm for psm in psms if (psm.hybrid_xcorr > psm.native_xcorr)]
        else:
            raise ValueError(f"psm_type must be '{NATIVE}' or '{HYBRID}'")

    @staticmethod
    def align_psms_to_proteome(
        psms: List["SpectrumPSMs"], fasta: Path
    ) -> List[HybridPosition]:
        fasta = Fasta(path=fasta)
        hybrid_positions = []
        for psm in psms:
            for hy in psm.hybrids:

                prot_pairs = list(
                    itertools.product(hy.left_proteins, hy.right_proteins)
                )
                for left_prot, right_prot in prot_pairs:
                    # Get where in protein it appears
                    left_pos = get_positions_of_subseq_in_seq(
                        subseq=hy.left_seq, seq=fasta.protein_name_to_seq_map[left_prot]
                    )
                    assert len(left_pos) == 1
                    left_pos = left_pos[0]
                    right_pos = get_positions_of_subseq_in_seq(
                        subseq=hy.right_seq,
                        seq=fasta.protein_name_to_seq_map[right_prot],
                    )
                    assert len(right_pos) == 1
                    right_pos = right_pos[0]
                    assert (
                        fasta.protein_name_to_seq_map[left_prot][
                            left_pos.inclusive_start : left_pos.exclusive_end
                        ]
                        == hy.left_seq
                    )
                    assert (
                        fasta.protein_name_to_seq_map[right_prot][
                            right_pos.inclusive_start : right_pos.exclusive_end
                        ]
                        == hy.right_seq
                    )

                    hybrid_positions.append(
                        HybridPosition(
                            left=ProteinRange.from_pos(protein=left_prot, pos=left_pos),
                            right=ProteinRange.from_pos(
                                protein=right_prot, pos=right_pos
                            ),
                        )
                    )
        return hybrid_positions


def get_pileup_from_positions(positions: List[ProteinRange]):
    pileup = defaultdict(lambda: defaultdict(int))
    for pos in positions:
        for idx in range(pos.inclusive_start, pos.exclusive_end):
            pileup[pos.protein][idx] += 1
    return dict(pileup)


def plot_hybrid_psm_pileup(
    psms: List[SpectrumPSMs], fasta: Path, out_path: Optional[Path] = None
):
    # Create pileups
    hybrid_positions = SpectrumPSMs.align_psms_to_proteome(psms=psms, fasta=fasta)

    left_pileup = get_pileup_from_positions(
        positions=[pos.left for pos in hybrid_positions]
    )
    right_pileup = get_pileup_from_positions(
        positions=[pos.right for pos in hybrid_positions]
    )
    prot_names = set(left_pileup.keys()).union(set(right_pileup.keys()))
    print(f"There are {len(prot_names)} proteins")

    # Plot the pileups for each protein
    fasta = Fasta(path=fasta)
    _, axs = fig_setup(nrows=len(prot_names), ncols=1, w=8)
    l, r = "left", "right"
    pileups = {l: left_pileup, r: right_pileup}
    colors = {l: "blue", r: "red"}
    max_cnt = max(
        max(flatten_list_of_lists(pileup.values() for pileup in left_pileup.values())),
        max(flatten_list_of_lists(pileup.values() for pileup in right_pileup.values())),
    )
    for idx, prot_name in enumerate(prot_names):
        ax = axs[idx]
        prot_seq = fasta.protein_name_to_seq_map[prot_name]
        for name, pileup in pileups.items():
            if prot_name not in pileup:
                continue
            data = np.array(
                [(ii, pileup[prot_name][ii]) for ii in range(len(prot_seq))]
            )
            data = data[data[:, 1] != 0]  # remove zeros
            _ = ax.plot(
                data[:, 0],
                data[:, 1],
                "o",
                color=colors[name],
                ms=2,
            )
        set_title_axes_labels(ax=ax, title=prot_name)
        _ = ax.set_ylim(bottom=0, top=max_cnt + 1)
        _ = ax.set_xlim(left=0, right=len(prot_seq) + 1)

    finalize(axs=axs)
    if out_path is not None:
        save_fig(
            path=out_path,
        )


@dataclass
class NeoFusionIteration:
    q: float
    min_score_delta: float
    fpr: float
    tp: int
    min_hybrid_score: float
    accepted_psm: List[str]


@dataclass
class NeoFusionRunner:
    q_vals: List[float]
    score_deltas: List[float]
    fpr_threshold: float

    @staticmethod
    def prepare_spectra_psms_for_neofusion(
        psms: List[SpectrumPSMs],
    ) -> pd.DataFrame:
        # Remove those PSMs that do not have a hybrid sequence or do not have a native q-value
        psms = SpectrumPSMs.filter_to_has_native_q(
            psms=SpectrumPSMs.filter_to_has_hybrid(psms=psms)
        )
        # Create dataframe for NeoFusion analysis
        df = pd.DataFrame(
            [
                [
                    psm.spectrum_uid,
                    psm.native_xcorr,
                    psm.hybrid_xcorr,
                    psm.native_q,
                ]
                for psm in psms
            ],
            columns=[
                "uid",
                "n_score",
                "h_score",
                "n_q",
            ],
        )
        df["delta"] = df.h_score - df.n_score
        return df

    @staticmethod
    def neofusion_iteration(
        df: pd.DataFrame,
        q_val: float,
        score_delta: float,
        fpr_thresh: float,
    ) -> Optional[NeoFusionIteration]:
        for colm in ["n_score", "h_score", "n_q", "delta"]:
            assert colm in df.columns
        neo_df = deepcopy(df)

        # Remove rows where hybrid score isn't high enough compared to native score
        neo_df = neo_df[neo_df.delta >= score_delta].copy()

        # Set which native PSMs are "gold-standard"
        neo_df["gold"] = neo_df.n_q <= q_val
        neo_df["fp"] = neo_df["gold"].copy()

        # Min hybrid score is lowest gold-standard native score
        min_hybrid_score = neo_df[neo_df.gold].n_score.min()
        if pd.isna(min_hybrid_score):
            return None

        # Sort in descending order
        neo_df.sort_values(by="h_score", ascending=False, inplace=True)
        neo_df.reset_index(drop=True, inplace=True)

        # Iterate through rows and, for each row i, find number number of false and true positives
        # in rows 1, 2, ..., i.
        fpr_colm = []
        tp_colm = []
        for row_idx, row in neo_df.iterrows():
            if row.h_score < min_hybrid_score:
                fpr_colm.append(None)
                tp_colm.append(None)
                continue

            tmp = neo_df.iloc[: row_idx + 1, :]
            fpr_colm.append(tmp.fp.sum() / tmp.shape[0])
            tp_colm.append((~tmp.fp).sum())
            # if not row.gold:
            #     # If hybrid score is greater than min hybrid score + score delta, set as gold
            #     if row.h_score >= (min_hybrid_score + score_delta):
            #         neo_df.at[row_idx, "gold"] = True
        neo_df["fpr"] = fpr_colm
        neo_df["tp"] = tp_colm
        if neo_df[neo_df["fpr"] < fpr_thresh].shape[0] > 0:
            # Get row that maximizes the number of true positives
            try:
                tmp = neo_df[neo_df["fpr"] < fpr_thresh]
                tp_maximizing_idx = tmp.tp.idxmax()
                tp_maximizing_row = tmp.loc[tp_maximizing_idx]
                accepted_psm = neo_df.uid.iloc[: tp_maximizing_idx + 1].tolist()

                return NeoFusionIteration(
                    q=q_val,
                    min_score_delta=score_delta,
                    fpr=tp_maximizing_row.fpr,
                    tp=tp_maximizing_row.tp,
                    min_hybrid_score=min_hybrid_score,
                    accepted_psm=accepted_psm,
                )
            except:
                print(f"Issue with q_val={q_val}, score_delta={score_delta}")
        return None

    def run_neofusion(
        self, psms: List[SpectrumPSMs]
    ) -> Tuple[NeoFusionIteration, List[SpectrumPSMs]]:
        df = self.prepare_spectra_psms_for_neofusion(psms=psms)

        neofusion_results = []
        for q_val in self.q_vals:
            for min_score_delta in self.score_deltas:
                result = self.neofusion_iteration(
                    df=df,
                    q_val=q_val,
                    score_delta=min_score_delta,
                    fpr_thresh=self.fpr_threshold,
                )
                if result is not None:
                    neofusion_results.append(result)

        best_iteration = max(neofusion_results, key=lambda x: x.tp)
        accepted_psm = list(
            filter(lambda psm: psm.spectrum_uid in best_iteration.accepted_psm, psms)
        )
        logger.info(
            f"Via NeoFusion, number of accepted hybrid PSM: {len(accepted_psm)}"
        )
        return best_iteration, accepted_psm


@dataclass
class AcceptedHybridPSMs:
    accepted_psms: List[SpectrumPSMs]

    @property
    def seq_to_accepted_psms(self) -> Dict[str, List[SpectrumPSMs]]:
        seq_to_psms = defaultdict(list)
        for psm in self.accepted_psms:
            seq_to_psms[psm.hybrid_seq].append(psm)
        return dict(seq_to_psms)

    @property
    def summary_dict(self):
        return {
            "num_accepted_hybrid_psms": len(self.accepted_psms),
            "num_accepted_unique_hybrid_seqs": len(self.seq_to_accepted_psms),
            "accepted_hybrid_psms": [psm.spectrum_uid for psm in self.accepted_psms],
        }


@dataclass
class TrueHybridSpectrumPSMsComparison:
    true_hybrid_seqs: Set[str]
    supported_true_seqs: List[Dict]

    @property
    def num_spectra_supporting_true_hybrids(self) -> int:
        return sum(
            seq_data["num_supporting_hybrid_psms"]
            for seq_data in self.supported_true_seqs
        )

    @classmethod
    def from_accepted_psms(
        cls, accepted_psms: AcceptedHybridPSMs, true_hybrid_seqs: Set[str]
    ) -> Dict:
        supported_true_seqs = []
        for seq in set(accepted_psms.seq_to_accepted_psms.keys()).intersection(
            true_hybrid_seqs
        ):
            seq_psms = accepted_psms.seq_to_accepted_psms[seq]
            data = {
                "seq": seq,
                "num_supporting_hybrid_psms": len(seq_psms),
                "spectra_uids": [psm.spectrum_uid for psm in seq_psms],
            }
            supported_true_seqs.append(data)
        return cls(
            true_hybrid_seqs=true_hybrid_seqs,
            supported_true_seqs=supported_true_seqs,
        )

    @property
    def summary_dict(self):
        return {
            "num_true_hybrid_seqs": len(self.true_hybrid_seqs),
            "num_supported_true_hybrid_seqs": len(self.supported_true_seqs),
            "num_spectra_supporting_true_hybrids": self.num_spectra_supporting_true_hybrids,
            "supported_true_hybrid_seqs": self.supported_true_seqs,
        }


@dataclass
class Experiment:
    name: str
    results_dir: Path
    hs_config: Path
    psms_path: Path

    def __post_init__(self):
        assert self.hs_config.exists(), f"HS config does not exist at {self.hs_config}"
        assert self.psms_path.exists(), f"PSMs path does not exist at {self.psms_path}"
        self.results_dir.mkdir(parents=True, exist_ok=True)

    @cached_property
    def psms(self):
        return SpectrumPSMs.load(path=self.psms_path)


def accept_hybrids_that_beat_natives(psms: List[SpectrumPSMs]) -> List[SpectrumPSMs]:
    accepted_psms = SpectrumPSMs.get_winners_by_xcorr(psms=psms, psm_type=HYBRID)
    logger.info(
        f"Number of accepted hybrid PSMs that beat natives by XCorr: {len(accepted_psms)}"
    )
    return accepted_psms


def accept_hybrid_psms_via_interpolated_q_value(
    psms: List[SpectrumPSMs], q_val_threshold: float = Q_VAL_THRESH
) -> List[SpectrumPSMs]:
    # Get those spectra that have a winning hybrid and the native has a q-value
    psms = SpectrumPSMs.filter_to_has_native_q(
        psms=SpectrumPSMs.get_winners_by_xcorr(psms=psms, psm_type=HYBRID)
    )
    native_q_interpolator = fit_native_xcorr_qval_interpolator(
        psms=SpectrumPSMs.filter_to_has_native_q(psms=psms)
    )
    accepted_psms = []
    for psm in psms:
        hybrid_q = native_q_interpolator(psm.hybrid_xcorr)
        if hybrid_q <= q_val_threshold:
            accepted_psms.append(psm)
    logger.info(
        f"Via interpolated q-value, number of accepted hybrid PSM: {len(accepted_psms)}"
    )
    return accepted_psms


def accept_hybrid_psms_via_neo_fusion(
    psms: List[SpectrumPSMs],
    q_vals: List[float] = DEFAULT_Q_RANGE,
    score_deltas: List[float] = DEFAULT_SCORE_CHANGE_RANGE,
    fpr_threshold: float = DEFAULT_FPR,
) -> List[SpectrumPSMs]:
    best_iteration, accepted_psms = NeoFusionRunner(
        q_vals=q_vals,
        score_deltas=score_deltas,
        fpr_threshold=fpr_threshold,
    ).run_neofusion(psms=psms)
    logger.info(f"Best NeoFusion iteration was: {best_iteration}")
    return accepted_psms


def psm_score_histogram(ax, data, label):
    _ = sns.histplot(
        data,
        # element="step",
        kde=True,
        stat="density",
        common_norm=False,
        ax=ax,
        label=f"{label} (n={len(data)})",
        # fill=False,
        alpha=0.4,
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


def remove_non_hybrid_psms(psms: List[CometPSM]) -> List[CometPSM]:
    logger.debug("Remove non-hybrid PSMs...")
    return [psm for psm in psms if psm.is_hybrid]


def get_protein_abundance_from_spectra_results(
    spectra_results: List[SpectrumCometResults],
    q_val_thresh: Literal[Q_VAL_THRESH] = Q_VAL_THRESH,
):
    assign_conf_psms = [
        spec_results.assign_conf
        for spec_results in spectra_results
        if spec_results.assign_conf is not None
    ]
    return ProteinAbundance.from_comet_psms(
        psms=assign_conf_psms,
        q_val_thresh=q_val_thresh,
    )


def add_protein_abundances(
    df: pd.DataFrame, hs_config: HypedsearchRunConfig, q_val_thresh: float
) -> pd.DataFrame:
    # Get protein abundances from native run
    native_results = SpectrumCometResults.from_dir(
        path=hs_config.native_run_dir, by_spectrum=False
    )
    prot_ab = ProteinAbundance.from_comet_psms(
        psms=[res.assign_conf for res in native_results if res.assign_conf is not None],
        q_val_thresh=q_val_thresh,
    )

    # For each row with a hybrid target, get the max abundance of the left and right proteins
    left_prot_abs = []
    left_prot_rel_abs = []
    right_prot_abs = []
    right_prot_rel_abs = []
    for row_idx, row in df.iterrows():
        if not pd.isna(row["hybrid_target_seq"]):
            left_prots = row["left_prots"].split(";")
            left_ab = max(prot_ab.get_ab(protein=prot) for prot in left_prots)
            left_rel_ab = max(prot_ab.get_rel_ab(protein=prot) for prot in left_prots)

            right_prots = row["right_prots"].split(";")
            right_ab = max(prot_ab.get_ab(protein=prot) for prot in right_prots)
            right_rel_ab = max(prot_ab.get_rel_ab(protein=prot) for prot in right_prots)
        else:
            left_ab = None
            right_ab = None
            left_rel_ab = None
            right_rel_ab = None
        left_prot_abs.append(left_ab)
        left_prot_rel_abs.append(left_rel_ab)
        right_prot_abs.append(right_ab)
        right_prot_rel_abs.append(right_rel_ab)
    df["left_prot_ab"] = left_prot_abs
    df["left_prot_rel_ab"] = left_prot_rel_abs
    df["right_prot_ab"] = right_prot_abs
    df["right_prot_rel_ab"] = right_prot_rel_abs
    return df


def plot_native_vs_hybrid_scores(
    psms: List[SpectrumPSMs],
    # true_hybrid_seqs: Set[str],
    score: str = XCORR,
    ax: Optional[Axes] = None,
) -> Axes:

    # Get SpectrumPSM objects with both a native and a hybrid target PSM
    native_and_hybrid_psms = [
        psm
        for psm in psms
        if (psm.native_target is not None) and (psm.hybrid_target is not None)
    ]

    if ax is None:
        _, axs = fig_setup()
        ax = axs[0]
    s = 7
    sns.scatterplot(
        x=[getattr(psm.native_target, score) for psm in native_and_hybrid_psms],
        y=[getattr(psm.hybrid_target, score) for psm in native_and_hybrid_psms],
        s=s,
        marker="o",
        color="blue",
        label=f"n={len(native_and_hybrid_psms)}",
        ax=ax,
    )
    # sns.scatterplot(
    #     data=data,
    #     x=f"native_target_{score}",
    #     y=f"hybrid_target_{score}",
    #     s=2 * s,
    #     marker="x",
    #     color="red",
    #     label=f"True hybrid seqs (n={data.shape[0]})",
    #     ax=ax,
    #     linewidth=1.5,
    # )
    plot_line(ax=ax, label="y=x")
    set_title_axes_labels(
        ax=ax,
        xlabel=f"Native target {score}",
        ylabel=f"Hybrid target\n{score}",
    )
    return ax


def fit_native_xcorr_qval_interpolator(
    psms: List[SpectrumPSMs],
    ax: Optional[Axes] = None,
) -> PchipInterpolator:
    # Get data
    df = pd.DataFrame(
        [(psm.native_xcorr, psm.native_q) for psm in psms if psm.native_q is not None],
        columns=[XCORR, Q_VAL],
    )
    xy = df.drop_duplicates(subset=XCORR)
    xy.sort_values(by=XCORR, inplace=True)
    x = xy[XCORR].to_numpy()
    y = xy[Q_VAL].to_numpy()
    interpolator = PchipInterpolator(x, y)

    # Plotting
    if ax is not None:
        # Plot data
        _ = sns.scatterplot(
            data=df, x=XCORR, y=Q_VAL, s=10, ax=ax, label=f"n={df.shape[0]}"
        )
        # Plot interpolating function
        x_new = np.linspace(x.min(), x.max(), 500)
        _ = ax.plot(x_new, interpolator(x_new), "r--", label="Interpolating function")
        ax.set_yscale("logit")
        set_title_axes_labels(
            ax=ax,
            xlabel=XCORR,
            ylabel=Q_VAL,
        )
        # On the right hand side axis, add the #(PSM with XCORR >= x)
        ax_copy = ax.twinx()
        y = [sum(df[XCORR] > x) for x in x_new]
        _ = sns.scatterplot(
            x=x_new,
            y=y,
            ax=ax_copy,
            # label="hybrid score",
            color="tab:red",
            s=7,
        )
        ax_copy.set_ylabel("Num PSMs with xcorr > x", color="tab:red")
        ax_copy.tick_params(axis="y", labelcolor="tab:red")
    return interpolator


def filter_to_top_n_highest_precursor_intensity_psms_per_mz(
    psms: List[SpectrumPSMs], n: int
) -> List[SpectrumPSMs]:
    mz_to_psms = defaultdict(list)
    for psm in psms:
        mz_to_psms[psm.spectrum.precursor_mz].append(psm)
    filtered_psms = []
    for mz, psms in mz_to_psms.items():
        psms = sorted(
            psms, key=lambda psm: psm.spectrum.precursor_abundance, reverse=True
        )
        filtered_psms.extend(psms[:n])
    return filtered_psms


def precursor_mz_plots(psms: List[SpectrumPSMs], sample: str, plot_dir: Path):
    mz_to_psms = defaultdict(list)
    for psm in psms:
        mz_to_psms[psm.spectrum.precursor_mz].append(psm)
    mz_to_psms = dict(mz_to_psms)

    num_supporting_spectra = Counter(len(psms) for mz, psms in mz_to_psms.items())
    _, axs = fig_setup(
        # nrows=1, ncols=2,
        w=10
    )
    ax = axs[0]
    _ = plot_sorted_1d_data(
        data=num_supporting_spectra,
        ax=ax,
        pt_labels=num_supporting_spectra,
        sort_idx=0,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="x = num supporting spectra",
        ylabel="Number of m/z's with num\nsupporting spectra = x",
        title=(
            f"sample = {sample}\n" + f"Num unique m/z values = {len(mz_to_psms)}\n"
            f"Num spectra = {sum(len(psms) for psms in mz_to_psms.values())}"
        ),
    )
    finalize(axs)
    save_fig(
        path=plot_dir / f"number_spectra_supporting_each_uniq_mz.png",
    )

    _, axs = fig_setup()
    cdf = ECDF(list(num_supporting_spectra.values()))
    ax = axs[0]
    _ = ax.plot(
        cdf.x,
        1 - cdf.y,
        "o-",
        ms=2,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="x = num supporting spectra",
        ylabel="Proportion of m/z's with\nnum supporting spectra >= x",
        title=(
            f"sample = {sample}\n" + f"Num unique m/z values = {len(mz_to_psms)}\n"
            f"num spectra = {sum(len(psms) for psms in mz_to_psms.values())}"
        ),
    )
    finalize(axs)
    save_fig(
        path=plot_dir / f"prop_of_uniq_mz_with_supporting_spectra.png",
    )


@click.command(
    name="spectrum-psms",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    Create SpectrumPSMs objects for Hypedsearch experiment
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    multiple=True,
    help="Paths to one or more Hypedsearch JSON configs",
)
@click.option(
    "--min_side_len",
    "-m",
    type=int,
    required=False,
    show_default=True,
    default=DEFAULT_MIN_SIDE_LEN,
    help="Minimum hybrid side length",
)
def cli_create_spectrum_psms(
    config: Tuple[Path, ...],
    min_side_len: int,
):
    for hs_config in config:
        # Creating SpectrumPSMs
        spectra_psms = SpectrumPSMs.from_hs_config(
            hs_config=hs_config, min_side_len=min_side_len, save=True
        )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_create_spectrum_psms)
    cli()
