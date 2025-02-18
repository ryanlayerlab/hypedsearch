import os
import re
import sqlite3
import subprocess
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Dict, Generator, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pydantic import BaseModel, Field, model_validator
from pyteomics import mzml

# Constants
## Directories and paths
from src.comet_utils import read_comet_txt_to_df
from src.erik_constants import (
    COMET_RUN_1_DIR,
    COMET_RUN_2_DIR,
    EXCLUSIVE_END,
    HS_DIR,
    INCLUSIVE_START,
    KMERS,
    MASS,
    MAX_KMER_LEN,
    PLAIN_PEPTIDE,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SAMPLE,
    SCAN,
    SEQ,
    SPECTRA_DIR,
    SPECTRUM_ID,
    THOMAS_SAMPLES,
)
from src.erik_utils import (
    compute_peptide_neutral_mass,
    compute_theoretical_mass_over_charge,
    dataclasses_to_df,
    pydantic_models_to_df,
)
from src.hypedsearch_utils import HSResult
from src.lookups.constants import (
    AMINO_ACID_MASSES,
    DOUBLY_CHARGED_B_BASE,
    DOUBLY_CHARGED_Y_BASE,
    PROTON_MASS,
    SINGLY_CHARGED_B_BASE,
    SINGLY_CHARGED_Y_BASE,
    WATER_MASS,
)
from src.lookups.data_classes import (
    Ion,
    Kmer,
    KmerWithMass,
    Peak,
    PeakIonMatch,
    Spectrum,
)
from src.lookups.protein_product_ion_db import (
    compute_b_ion_neutral_mass,
    compute_y_ion_neutral_mass,
)
from src.plot_utils import fig_setup, finalize, set_title_axes_labels


def parse_mzml(mzml_path: Union[str, Path]) -> Generator[Spectrum, None, None]:
    mzml_path = Path(mzml_path)
    spectra = mzml.read(str(mzml_path))
    for spectrum in spectra:
        yield parse_spectrum(spectrum=spectrum, mzml=mzml_path.stem)


def parse_spectrum(spectrum: Dict, mzml: Optional[str] = None) -> Spectrum:
    masses, abundances = tuple(spectrum["m/z array"]), tuple(
        spectrum["intensity array"]
    )
    peaks = [
        Peak(mz=masses[idx], abundance=abundances[idx], id=idx)
        for idx in range(len(masses))
    ]
    spectrum_id = spectrum.get("id", "")
    scan_num = int(spectrum_id.split("=")[1])
    return Spectrum(
        # mass_over_charges=masses,
        # abundances=abundances,
        scan_num=scan_num,
        peaks=peaks,
        precursor_mz=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
            "selectedIon"
        ][0]["selected ion m/z"],
        precursor_charge=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
            "selectedIon"
        ][0]["charge state"],
        precursor_abundance=spectrum["precursorList"]["precursor"][0][
            "selectedIonList"
        ]["selectedIon"][0]["peak intensity"],
        spectrum_id=spectrum_id,
        retention_time=spectrum["scanList"]["scan"][0]["scan start time"],
        mzml=mzml,
    )


def plot_mass_spectrum(spectrum: Spectrum):

    _, axs = fig_setup(nrows=1, ncols=1)
    mz_values = list(spectrum.mass_over_charges)
    abundances = list(spectrum.abundances)

    # _ = axs[0].stem(mz_values, abundances)
    _ = axs[0].bar(mz_values, abundances, width=5)
    _ = set_title_axes_labels(axs[0], xlabel="m/z", ylabel="intensity")
    finalize(axs)
    # plt.tight_layout()
    return axs


def relative_abundance_filtering(mass_over_charges, abundances, rel_abundance_cutoff):
    """
    Selects only those ions whose relative abundance is above a given percentage
    of the highest abundance
    Relative abundance measures the intensity of a given ion signal relative to the
    intensity of the most abundant ion.
    """
    assert (rel_abundance_cutoff >= 0) and (
        rel_abundance_cutoff <= 1
    ), f"The relative abundance cutoff proportion should be between 0 and 1. Received {rel_abundance_cutoff}"
    relative_val = max(abundances)
    min_value = relative_val * rel_abundance_cutoff

    mass_abundances = zip(mass_over_charges, abundances)
    filtered_mass_abundances = [x for x in mass_abundances if x[1] >= min_value]
    mass_over_charges = tuple(float(x) for x, _ in filtered_mass_abundances)
    abundances = tuple(float(x) for _, x in filtered_mass_abundances)

    return (mass_over_charges, abundances)


def get_mass_matching_kmers(
    mass: float,
    db_cursor,
    tolerance: float = 10,
    query_type: Literal["count", "both", "seq", "mass"] = "count",
):
    query = get_query_to_select_rows_by_mass(
        mass=mass, mass_tol=tolerance, query_type=query_type
    )
    db_cursor.execute(query)
    return db_cursor.fetchall()


def get_query_to_select_rows_by_mass(
    mass: float,
    mass_tol: float,
    query_type: Literal["count", "both", "seq", "mass"] = "count",
):
    if query_type == "count":
        query = f"SELECT COUNT(*) from {KMERS} WHERE {MASS} BETWEEN {mass - mass_tol} AND {mass + mass_tol}"
    elif query_type == "both":
        query = f"SELECT {SEQ}, {MASS} FROM kmers WHERE {MASS} BETWEEN {mass - mass_tol} AND {mass + mass_tol}"
    elif query_type == "seq":
        query = f"SELECT {SEQ} FROM kmers WHERE {MASS} BETWEEN {mass - mass_tol} AND {mass + mass_tol}"
    elif query_type == "mass":
        query = f"SELECT {MASS} FROM kmers WHERE {MASS} BETWEEN {mass - mass_tol} AND {mass + mass_tol}"
    return query


def query_database(query, db_path):
    # Open a new connection for each thread
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    conn.close()
    return result


# Process peaks
# if filter_num_peaks is not None:
#         masses, abundances = select_highest_abundance_peaks(
#             mass_to_charges=masses, abundances=abundances
#         )
#     if relative_abundance_filter_percentage is not None:
#         masses, abundances = relative_abundance_filtering(
#             mass_to_charges=masses,
#             abundances=abundances,
#             rel_abundance_cutoff=relative_abundance_filter_percentage / 100,
#         )
@dataclass
class CombinedData:
    spectrum: pd.DataFrame
    comet_1: pd.DataFrame
    comet_2: pd.DataFrame
    hs: pd.DataFrame


def spectra_to_dataframe(spectra: List[Spectrum]) -> pd.DataFrame:
    return pd.DataFrame([asdict(spectrum) for spectrum in spectra])


def get_specific_spectrum(sample: str, scan_num: int) -> Spectrum:
    """
    Args:
        - scan_num is the spectrum's 1-based index
    """
    mzml_path = SPECTRA_DIR / f"{sample}.mzML"
    matched_spectrum = None
    spectra = parse_mzml(mzml_path)
    matched_spectrum = None
    for spectrum in spectra:
        if spectrum.scan_num == scan_num:
            matched_spectrum = spectrum
            break
    return matched_spectrum


def get_data_for_spectrum(sample: str, scan_num: int):
    """
    Args:
        - scan_num is the spectrum's 1-based index
    """
    spectrum_idx = scan_num - 1

    # Get spectrum from MZML
    spectrum = get_specific_spectrum(sample=sample, scan_num=scan_num)

    # Get Comet run 1 data
    comet_output = COMET_RUN_1_DIR / f"{sample}/{sample}.txt"
    comet_1_df = read_comet_txt_to_df(txt_path=comet_output)

    # Get Comet run 2 data
    comet_output = COMET_RUN_2_DIR / f"{sample}.txt"
    comet_2_df = read_comet_txt_to_df(txt_path=comet_output)

    # Get hypedsearch data
    hs_path = HS_DIR / f"HS_{sample}.txt"
    hs_df = pydantic_models_to_df(HSResult.from_txt(file_path=hs_path))

    return CombinedData(
        spectrum=spectrum,
        comet_1=comet_1_df[comet_1_df[SCAN] == scan_num],
        comet_2=comet_2_df[comet_2_df[SCAN] == scan_num],
        hs=hs_df[hs_df[SPECTRUM_ID] == spectrum_idx],
    )


def get_b_ion_sequences(peptide: str) -> List[str]:
    return [peptide[:i] for i in range(1, len(peptide) + 1)]  # Prefixes (b-ions)


def get_y_ion_sequences(peptide: str) -> List[str]:
    return [peptide[i:] for i in range(1, len(peptide))]  # Suffixes (y-ions)


def get_b_y_ion_sequences(peptide: str):
    b_ions = get_b_ion_sequences(peptide=peptide)
    y_ions = get_y_ion_sequences(peptide=peptide)
    return b_ions + y_ions


def generate_kmers(peptide: str, max_k: int, min_k: Optional[int] = 1) -> List[Kmer]:
    """"""
    disallowed_amino_acid_symbols = ["B", "X", "U", "Z", "O", "J"]
    kmers = []
    for kmer_len in range(min_k, max_k + 1):  # Iterate over all k from 1 to max_k
        for start in range(len(peptide) - kmer_len + 1):  # Slide over the string
            end = start + kmer_len
            kmer = peptide[start:end]
            if any(char in disallowed_amino_acid_symbols for char in kmer):
                continue
            kmers.append(
                Kmer(
                    seq=peptide[start : start + kmer_len],
                    inclusive_start=start,
                    exclusive_end=end,
                )
            )

    return kmers


#     assert charge <= 2, f"Only charge=1 or charge=2 are supported. Got charge={charge}"
#     if ion == "y":
#         total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
#         total += sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
#     elif ion == "b":
#         total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
#         total += sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
#     else:
#         raise ValueError(f"Ion must be 'b' or 'y'. Got ion={ion}")
#     mz = total / charge
#     return mz


def generate_kmers_with_masses(
    peptide: str,
    max_kmer_len: int = MAX_KMER_LEN,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    charge: int = 1,
    # peptide_id: Optional[int] = None,
) -> Generator[Kmer, None, None]:
    kmers = generate_kmers(peptide=peptide, max_k=max_kmer_len)
    for kmer in kmers:
        yield KmerWithMass(
            kmer=kmer,
            mass=compute_peptide_neutral_mass(
                peptide=kmer.seq,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
                charge=charge,
            ),
        )


def get_indices_of_largest_elements(array: List[float], top_n: int):
    if top_n >= len(array):
        return np.arange(0, len(array))
    array = np.array(array)
    # Get the indices of the largest N elements
    indices = np.argpartition(-array, top_n)[:top_n]
    # Sort these indices to have them in descending order of the elements
    sorted_indices = np.sort(indices)
    return sorted_indices


def mass_difference_in_ppm(ref_mass: float, query_mass: float) -> float:
    return (abs(ref_mass - query_mass) / ref_mass) * (10**6)


def load_comet_data(samples: List[str] = THOMAS_SAMPLES, run: int = 1):
    if run == 1:
        comet_results_dir = COMET_RUN_1_DIR
    elif run == 2:
        comet_results_dir = COMET_RUN_2_DIR
    else:
        raise ValueError(f"'run' must be 1 or 2. You set run={run}")
    comet_dfs = []
    for sample in samples:
        print(f"Reading data for {sample}")
        comet_output = comet_results_dir / f"{sample}/{sample}.txt"
        comet_dfs.append(read_comet_txt_to_df(txt_path=comet_output))
    comet_df = pd.concat(comet_dfs, ignore_index=True)
    return comet_df


def load_mzml_data(samples: List[str] = THOMAS_SAMPLES):
    mzml_data = []
    for sample in samples:
        print(f"Reading sample {sample}'s MZML")
        mzml_path = SPECTRA_DIR / f"{sample}.mzML"
        spectra = parse_mzml(mzml_path=mzml_path)
        mzml_data.extend(list(spectra))
    # if is_dataclass(mzml_data[0]):
    #     df = dataclasses_to_df(dataclasses=mzml_data)
    # else:
    #     df = pydantic_models_to_df(models=mzml_data)

    # # Add scan number column by extracting the number in "scan=<N>"
    # df[SCAN] = df["spectrum_id"].apply(
    #     lambda id: int(re.search(r"scan=(\d+)", id).group(1))
    # )
    # return df
    return mzml_data


def calculate_expected_precursor_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
):
    mass = sum([amino_acid_mass_lookup[aa] for aa in aa_seq]) + WATER_MASS
    return (mass + (charge * PROTON_MASS)) / charge


@dataclass(frozen=True)
class PeakMatchingIons:
    b_matches: Optional[List[PeakIonMatch]] = None
    y_matches: Optional[List[PeakIonMatch]] = None


def top_n_peak_filtering(peaks: List[Peak], n: int) -> List[Peak]:
    abundances = [peak.abundance for peak in peaks]
    indices = get_indices_of_largest_elements(array=abundances, top_n=n)
    peaks = np.array(peaks)
    return list(peaks[indices])


def find_peaks_that_match_aa_seq(
    aa_seq: str, peaks: List[Peak], ppm_tolerance: float, charges_to_consider
) -> List[Peak]:
    matches = []
    total_num_ions = 0
    for _, peak in enumerate(peaks):
        for charge in charges_to_consider:
            # b-ions
            ion_mz = compute_b_ion_neutral_mass(aa_seq=aa_seq, charge=charge)
            total_num_ions
            ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mz)
            if ppm_diff <= ppm_tolerance:
                matches.append(
                    PeakIonMatch(
                        peak=peak,
                        ion=Ion(seq=aa_seq, charge=charge, mz=ion_mz, ion_type="b"),
                    )
                )

            # y-ions
            ion_mz = compute_y_ion_neutral_mass(aa_seq=aa_seq, charge=charge)
            ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mz)
            if ppm_diff <= ppm_tolerance:
                matches.append(
                    PeakIonMatch(
                        peak=peak,
                        ion=Ion(seq=aa_seq, charge=charge, mz=ion_mz, ion_type="y"),
                    )
                )
    return matches


def find_ions_that_match_peak(
    peptide: str, charges_to_consider: List[int], ppm_tolerance: float, peak: Peak
) -> List[Ion]:
    y_seqs = get_y_ion_sequences(peptide=peptide)
    b_seqs = get_b_ion_sequences(peptide=peptide)

    matches = []
    # b-ions
    for ion in b_seqs:
        for charge in charges_to_consider:
            ion_mz = compute_b_ion_neutral_mass(aa_seq=ion, charge=charge)
            ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mz)
            if ppm_diff <= ppm_tolerance:
                matches.append(Ion(seq=ion, charge=charge, mz=ion_mz, ion_type="b"))

    for ion in y_seqs:
        for charge in charges_to_consider:
            ion_mz = compute_y_ion_neutral_mass(aa_seq=ion, charge=charge)
            ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mz)
            if ppm_diff <= ppm_tolerance:
                matches.append(Ion(seq=ion, charge=charge, mz=ion_mz, ion_type="y"))

    return matches


def get_ion_matches_between_peptide_and_spectrum_peaks(
    peptide: str,
    charges_to_consider: List[int],
    ppm_tol: int,
    peaks: List[Peak],
) -> PeakMatchingIons:
    y_seqs = get_y_ion_sequences(peptide=peptide)
    b_seqs = get_b_ion_sequences(peptide=peptide)

    # b-ions
    b_matches = []
    for _, peak in enumerate(peaks):
        for ion in b_seqs:
            for charge in charges_to_consider:
                ion_mz = compute_b_ion_neutral_mass(aa_seq=ion, charge=charge)
                ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mz)
                if ppm_diff <= ppm_tol:
                    b_matches.append(
                        PeakIonMatch(
                            peak=peak,
                            ion=Ion(seq=ion, charge=charge, mz=ion_mz, ion_type="b"),
                        )
                    )

    # y-ions
    y_matches = []
    for _, peak in enumerate(peaks):
        for ion in y_seqs:
            for charge in charges_to_consider:
                ion_mass = compute_y_ion_neutral_mass(aa_seq=ion, charge=charge)
                ppm_diff = mass_difference_in_ppm(ref_mass=peak.mz, query_mass=ion_mass)
                if ppm_diff <= ppm_tol:
                    y_matches.append(
                        PeakIonMatch(
                            peak=peak,
                            ion=Ion(seq=ion, charge=charge, mz=ion_mz, ion_type="y"),
                        )
                    )

    return PeakMatchingIons(b_matches=b_matches, y_matches=y_matches)


def compare_comet_output_row_to_corresponding_spectrum(
    row: pd.Series,
    charges_to_consider: List[int],
    ppm_tolerance: float,
    top_n_peaks: Optional[int] = None,
):
    # Get spectrum
    spectrum = get_specific_spectrum(sample=row[SAMPLE], scan_num=row[SCAN])

    # Get ion matches
    matches = get_ion_matches_between_peptide_and_spectrum_peaks(
        spectrum=spectrum,
        peptide=row[PLAIN_PEPTIDE],
        charges_to_consider=charges_to_consider,
        ppm_tol=ppm_tolerance,
        top_n_peaks=top_n_peaks,
    )

    # Analyze matches
    peak_ids = []
    b_matches, y_matches = 0, 0
    for match in matches:
        peak_ids.append(match.peak.id)
        if match.ion.ion_type == "b":
            b_matches += 1
        else:
            y_matches += 1

    return {
        "num_matching_peaks": len(set(peak_ids)),
        "num_b_matches": b_matches,
        "num_y_matches": y_matches,
    }


@dataclass(frozen=True)
class PeptideSpectrumMatchingInfo:
    aa_seq: str
    total_num_peaks: int
    num_peaks_with_ion_match: int
    total_ion_seqs: int
    num_ion_seqs_with_peak_match: int
    charges_to_consider: List[int]
    top_n_peaks: Optional[int] = None


def compare_peptide_to_spectrum(
    spectrum: Spectrum,
    peptide: str,
    charges_to_consider: List[int],
    ppm_tolerance: float,
    top_n_peaks: Optional[int] = None,
):
    if top_n_peaks is not None:
        peaks = top_n_peak_filtering(peaks=spectrum.peaks, n=top_n_peaks)
    else:
        peaks = spectrum.peaks

    # Get information for each peak
    num_peaks_with_ion_match = 0
    for peak in peaks:
        matching_ions = find_ions_that_match_peak(
            peptide=peptide,
            charges_to_consider=charges_to_consider,
            ppm_tolerance=ppm_tolerance,
            peak=peak,
        )
        if len(matching_ions) > 0:
            num_peaks_with_ion_match += 1
        # print(peak)
        # print(matching_ions)

    # Get information for each b- and y-ion AA sequence
    ion_seqs = get_b_ion_sequences(peptide=peptide) + get_y_ion_sequences(
        peptide=peptide
    )
    num_ion_seqs_with_peak_match = 0
    for aa_seq in ion_seqs:
        peak_ion_matches = find_peaks_that_match_aa_seq(
            aa_seq=aa_seq,
            peaks=peaks,
            ppm_tolerance=ppm_tolerance,
            charges_to_consider=charges_to_consider,
        )
        if len(peak_ion_matches) > 0:
            num_ion_seqs_with_peak_match += 1

    return PeptideSpectrumMatchingInfo(
        aa_seq=peptide,
        total_num_peaks=len(peaks),
        num_peaks_with_ion_match=num_peaks_with_ion_match,
        total_ion_seqs=len(ion_seqs),
        num_ion_seqs_with_peak_match=num_ion_seqs_with_peak_match,
        charges_to_consider=charges_to_consider,
        top_n_peaks=top_n_peaks,
    )
