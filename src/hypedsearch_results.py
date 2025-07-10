from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path
from time import time
from typing import Counter, Dict, List, Literal, Optional, Set, Union

import numpy as np
import pandas as pd

from src.comet_utils import CometPSM, read_comet_txts_in_dir
from src.constants import HS_PREFIX
from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import Spectrum
from src.peptide_spectrum_comparison import PeptideSpectrumComparison
from src.utils import (
    dataclass_list_to_df,
    decompress_and_depickle,
    flatten_list_of_lists,
    get_rank,
    get_time_in_diff_units,
    pickle_and_compress,
)


@dataclass
class NativePSM(CometPSM):
    """
    Class to hold native PSMs
    """

    num_seqs: int
    num_seqs_supported: int
    prop_intensity_supported: float
    num_b_seqs: int
    num_b_seqs_supported: int
    b_ion_seqs_supported: Set[str]
    prop_intensity_supported_by_b_seqs: float
    num_y_seqs: int
    num_y_seqs_supported: int
    y_ion_seqs_supported: Set[str]
    prop_intensity_supported_by_y_seqs: float
    # seqs_df: pd.DataFrame

    @staticmethod
    def _create_native_creation_dict(comet_psm: CometPSM, spectrum: Spectrum) -> Dict:
        comp = PeptideSpectrumComparison.from_str(
            spectrum=spectrum, peptide=comet_psm.seq
        )
        tmp = asdict(comet_psm)
        tmp.update(
            {
                "num_seqs": comp.num_b_ion_seqs + comp.num_y_ion_seqs,
                "num_seqs_supported": comp.num_ion_seqs_supported,
                "prop_intensity_supported": comp.prop_of_intensity_supported_by_seqs,
                "num_b_seqs": comp.num_b_ion_seqs,
                "num_b_seqs_supported": comp.num_b_ion_seqs_supported,
                "b_ion_seqs_supported": list(comp.b_ion_seqs_supported),
                "prop_intensity_supported_by_b_seqs": comp.prop_of_intensity_supported_by_b_seqs,
                "num_y_seqs": comp.num_y_ion_seqs,
                "num_y_seqs_supported": comp.num_y_ion_seqs_supported,
                "y_ion_seqs_supported": list(comp.y_ion_seqs_supported),
                "prop_intensity_supported_by_y_seqs": comp.prop_of_intensity_supported_by_y_seqs,
                # "seqs_df": comp.product_ion_seqs_with_matching_peaks,
            }
        )
        return tmp

    @classmethod
    def from_native_comet_psm_and_spectrum(
        cls, comet_psm: CometPSM, spectrum: Spectrum
    ) -> "NativePSM":
        creation_dict = cls._create_native_creation_dict(
            comet_psm=comet_psm, spectrum=spectrum
        )
        return cls(**creation_dict)


def get_max_protein_count(proteins: List[str], protein_counts: Union[Path, Counter]):
    if isinstance(protein_counts, Path):
        protein_counts = decompress_and_depickle(protein_counts)
    return max([protein_counts[prot] for prot in proteins])


def get_protein_count_rank(proteins: List[str], protein_counts: Union[Path, Counter]):
    if isinstance(protein_counts, Path):
        protein_counts = decompress_and_depickle(protein_counts)
    max_count = get_max_protein_count(proteins=proteins, protein_counts=protein_counts)
    return get_rank(values=protein_counts.values(), query_val=max_count)


@dataclass
class HybridPSM(NativePSM):
    b_seq: str
    y_seq: str
    b_proteins: Set[str]
    y_proteins: Set[str]
    num_hybrid_seqs: int
    num_hybrid_seqs_supported: int
    hybrid_seqs_supported: Set[str]
    prop_intensity_supported_by_hybrid_seqs: float

    @classmethod
    def from_hybrid_comet_psm_and_spectrum(
        cls,
        comet_psm: CometPSM,
        spectrum: Spectrum,
        hybrid_seqs: Dict[str, List[HybridPeptide]],
    ) -> List["HybridPSM"]:
        native_dict = cls._create_native_creation_dict(
            comet_psm=comet_psm, spectrum=spectrum
        )
        hybrid_psms = []
        for hybrid_peptide in hybrid_seqs[comet_psm.seq]:
            comp = PeptideSpectrumComparison.from_str(
                spectrum=spectrum, peptide=hybrid_peptide.seq_with_hyphen
            )
            tmp = deepcopy(native_dict)
            tmp.update(
                {
                    "b_seq": hybrid_peptide.b_seq,
                    "y_seq": hybrid_peptide.y_seq,
                    "b_proteins": hybrid_peptide.b_prot_names,
                    "y_proteins": hybrid_peptide.y_prot_names,
                    "num_hybrid_seqs": comp.num_hybrid_seqs,
                    "num_hybrid_seqs_supported": comp.num_hybrid_seqs_supported,
                    "hybrid_seqs_supported": list(comp.hybrid_seqs_supported),
                    "prop_intensity_supported_by_hybrid_seqs": comp.prop_of_intensity_supported_by_hybrid_seqs,
                }
            )
            hybrid_psms.append(cls(**tmp))
        return hybrid_psms

    @staticmethod
    def get_right_side_count(
        y_seq: str,
        comet_psms: List[CometPSM],
    ) -> int:
        # Get Comet PSMs from which to calculate abundance
        right_side_matching_psms = [
            psm for psm in comet_psms if psm.seq.startswith(y_seq[: len(psm.seq)])
        ]
        return len(right_side_matching_psms)


@dataclass
class HSResult:
    """
    Class to hold the result of a HypedSearch run
    """

    spectrum: Spectrum
    native_psms: List[CometPSM]
    hybrid_psms: List[CometPSM]
    seq_to_hybrid_peptides: Dict[str, List[HybridPeptide]]
    hybrid_forming_info: Dict


def process_hypedsearch_results(
    hs_result: HSResult,
    output_type: Literal["csv", "pklz"],
    output_dir: Path,
):
    output_paths = hypedsearch_output_file_names(
        output_dir=output_dir,
        output_type=output_type,
        mzml_name=hs_result.spectrum.mzml.stem,
        scan=hs_result.spectrum.scan,
    )
    if output_type == "pklz":
        output_path = output_paths[0]  # there should just be one path
        pickle_and_compress(
            obj=hs_result,
            file_path=output_path,
        )
    elif output_type == "csv":
        # Save native PSMs
        native_psms = [
            NativePSM.from_native_comet_psm_and_spectrum(
                comet_psm=psm, spectrum=hs_result.spectrum
            )
            for psm in hs_result.native_psms
        ]
        if len(native_psms) > 0:
            native_df = dataclass_list_to_df(dataclass_list=native_psms)
            output_path = [path for path in output_paths if "native" in path.name][
                0
            ]  # should only be one
            native_df.to_csv(output_path, index=False)

        # Save hybrid PSMs
        hybrid_psms = flatten_list_of_lists(
            [
                HybridPSM.from_hybrid_comet_psm_and_spectrum(
                    comet_psm=psm,
                    spectrum=hs_result.spectrum,
                    hybrid_seqs=hs_result.seq_to_hybrid_peptides,
                )
                for psm in hs_result.hybrid_psms
            ]
        )
        if len(hybrid_psms) > 0:
            hybrid_df = dataclass_list_to_df(dataclass_list=hybrid_psms)
            output_path = [path for path in output_paths if "hybrid" in path.name][0]
            hybrid_df.to_csv(output_path, index=False)

    else:
        raise ValueError(
            "output_type must be either 'csv' or 'pklz'. " f"Got {output_type} instead."
        )


def hypedsearch_output_file_names(
    output_dir: Path, output_type: Literal["csv", "pklz"], mzml_name: str, scan: int
) -> List[Path]:
    if output_type == "pklz":
        paths = [output_dir / f"{HS_PREFIX}{mzml_name}_{scan}.pklz"]
    elif output_type == "csv":
        paths = [
            output_dir / f"{HS_PREFIX}{mzml_name}_{scan}_natives.csv",
            output_dir / f"{HS_PREFIX}{mzml_name}_{scan}_hybrids.csv",
        ]
    return paths


def add_counts(
    hs_results: Union[List[HSResult], List[Path]],
    comet_results_dir: Path,
    protein_counts: Path,
):
    comet_psms = read_comet_txts_in_dir(comet_results_dir)
    comet_psms = [psm for psm in comet_psms if psm.num == 1]  # only use top hits
    protein_counts = decompress_and_depickle(protein_counts)

    native_dfs, hybrid_dfs = [], []
    num_results = len(hs_results)
    run_times = []
    for idx, hs_result in enumerate(hs_results):
        start_time = time()
        if idx % 100 == 0:
            print(f"Processing {idx + 1}/{num_results} result")
            avg_run_time = np.mean(run_times)
            remaining_results = num_results - idx
            print(f"Avg run time: {get_time_in_diff_units(avg_run_time)} seconds.")
            print(f"ETA = {get_time_in_diff_units(remaining_results * avg_run_time)}")
        if isinstance(hs_result, Path):
            hs_result = decompress_and_depickle(hs_result)
        # Save native PSMs
        native_psms = [
            NativePSM.from_native_comet_psm_and_spectrum(
                comet_psm=psm, spectrum=hs_result.spectrum
            )
            for psm in hs_result.native_psms
        ]
        native_df = dataclass_list_to_df(dataclass_list=native_psms)
        native_df["max_protein_count"] = native_df.proteins.apply(
            lambda proteins: get_max_protein_count(
                proteins=proteins, protein_counts=protein_counts
            )
        )
        native_df["max_protein_rank"] = native_df.proteins.apply(
            lambda proteins: get_protein_count_rank(
                proteins=proteins, protein_counts=protein_counts
            )
        )
        native_dfs.append(native_df)

        # Save hybrid PSMs
        hybrid_psms = flatten_list_of_lists(
            [
                HybridPSM.from_hybrid_comet_psm_and_spectrum(
                    comet_psm=psm,
                    spectrum=hs_result.spectrum,
                    hybrid_seqs=hs_result.seq_to_hybrid_peptides,
                )
                for psm in hs_result.hybrid_psms
            ]
        )
        hybrid_df = dataclass_list_to_df(dataclass_list=hybrid_psms)
        hybrid_df["max_b_protein_count"] = hybrid_df.b_proteins.apply(
            lambda proteins: get_max_protein_count(
                proteins=proteins, protein_counts=protein_counts
            )
        )
        hybrid_df["max_y_protein_count"] = hybrid_df.y_proteins.apply(
            lambda proteins: get_max_protein_count(
                proteins=proteins, protein_counts=protein_counts
            )
        )
        hybrid_df["right_side_count"] = hybrid_df.y_seq.apply(
            lambda y_seq: HybridPSM.get_right_side_count(
                y_seq=y_seq,
                comet_psms=comet_psms,
            )
        )
        hybrid_dfs.append(hybrid_df)
        duration = time() - start_time
        run_times.append(duration)
    native_df = pd.concat(native_dfs, ignore_index=True)
    hybrid_df = pd.concat(hybrid_dfs, ignore_index=True)
    return native_df, hybrid_df
