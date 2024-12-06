import os
import sqlite3
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Literal, Optional, Tuple, Union

from matplotlib import pyplot as plt
from pydantic import BaseModel, Field, model_validator
from pyteomics import mzml

# Constants
## Directories and paths
from src.comet_utils import read_comet_txt_to_df
from src.erik_constants import *
from src.erik_utils import *
from src.hypedsearch_utils import HSResult
from src.lookups.constants import AMINO_ACID_MASSES
from src.plot_utils import fig_setup, finalize, set_title_axes_labels


# New classes
class Peak(BaseModel):
    mz: float
    abundance: float


class Spectrum(BaseModel):
    mass_over_charges: Tuple[float, ...]
    abundances: Tuple[float, ...]
    precursor_mass: float
    precursor_charge: int
    precursor_abundance: float
    spectrum_id: str
    retention_time: float
    peaks: Optional[Tuple[Peak]] = Field(init=False, default=None)

    def model_post_init(self, __context: dict) -> None:
        assert len(self.mass_over_charges) == len(self.abundances)
        self.peaks = tuple(
            (
                Peak(mz=self.mass_over_charges[ii], abundance=self.abundances[ii])
                for ii in range(len(self.mass_over_charges))
            )
        )


def parse_mzml(mzml_path: Union[str, Path]):
    spectra = mzml.read(str(mzml_path))
    for spectrum in spectra:
        yield parse_spectrum(spectrum=spectrum)


def parse_spectrum(
    spectrum: Dict,
):
    masses, abundances = tuple(spectrum["m/z array"]), tuple(
        spectrum["intensity array"]
    )
    return Spectrum(
        mass_over_charges=masses,
        abundances=abundances,
        precursor_mass=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
            "selectedIon"
        ][0]["selected ion m/z"],
        precursor_charge=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
            "selectedIon"
        ][0]["charge state"],
        precursor_abundance=spectrum["precursorList"]["precursor"][0][
            "selectedIonList"
        ]["selectedIon"][0]["peak intensity"],
        spectrum_id=spectrum.get("id", ""),
        retention_time=spectrum["scanList"]["scan"][0]["scan start time"],
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
        mass=mass, tolerance=tolerance, query_type=query_type
    )
    db_cursor.execute(query)
    return db_cursor.fetchall()


def get_query_to_select_rows_by_mass(
    mass: float,
    tolerance: float,
    query_type: Literal["count", "both", "seq", "mass"] = "count",
):
    if query_type == "count":
        query = f"SELECT COUNT(*) from {KMERS} WHERE {MASS} BETWEEN {mass - tolerance} AND {mass + tolerance}"
    elif query_type == "both":
        query = f"SELECT {SEQ}, {MASS} FROM kmers WHERE {MASS} BETWEEN {mass - tolerance} AND {mass + tolerance}"
    elif query_type == "seq":
        query = f"SELECT {SEQ} FROM kmers WHERE {MASS} BETWEEN {mass - tolerance} AND {mass + tolerance}"
    elif query_type == "mass":
        query = f"SELECT {MASS} FROM kmers WHERE {MASS} BETWEEN {mass - tolerance} AND {mass + tolerance}"
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


def get_specific_spectrum(sample: str, scan_num: int) -> pd.Series:
    """
    Args:
        - scan_num is the spectrum's 1-based index
    """
    mzml_path = SPECTRA_DIR / f"{sample}.mzML"
    spectra = pydantic_models_to_df(list(parse_mzml(mzml_path)))
    spectrum_of_interest = spectra[spectra[SPECTRUM_ID] == f"scan={scan_num}"]
    return spectrum_of_interest.iloc[0]


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


@dataclass
class Peptide:
    peptide: str

    @property
    def b_y_ions(self):
        return get_b_y_ions(peptide=self.peptide)

    @property
    def theoretical_mass(self, charge: int = 1):
        return compute_theoretical_mass_over_charge(peptide=self.peptide, charge=charge)

    @property
    def theoretical_b_y_ion_mass_spectrum(
        self, charge: int = 1, const_abundance: float = 1
    ):
        peaks = []
        for ion in self.b_y_ions:
            peaks.append(
                Peak(
                    mz=compute_theoretical_mass_over_charge(peptide=ion, charge=charge),
                    abundance=const_abundance,
                )
            )
        return peaks


def get_b_y_ions(peptide: str):
    b_ions = [peptide[:i] for i in range(1, len(peptide) + 1)]  # Prefixes (b-ions)
    y_ions = [peptide[i:] for i in range(1, len(peptide))]  # Suffixes (y-ions)
    return b_ions + y_ions


def generate_kmers(peptide: str, max_k: int):
    disallowed_amino_acid_symbols = ["B", "X", "U", "Z", "O", "J"]
    kmers = []
    for k in range(1, max_k + 1):  # Iterate over all k from 1 to max_k
        for start in range(len(peptide) - k + 1):  # Slide over the string
            kmer = peptide[start : start + k]
            if any(char in disallowed_amino_acid_symbols for char in kmer):
                continue
            kmers.append(peptide[start : start + k])  # Extract the k-mer
    return kmers


def generate_kmers_with_masses(
    peptide: str,
    max_kmer_len: int = MAX_PEPTIDE_LEN,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    charge: int = 1,
) -> Generator[Kmer]:
    """
    Assumes charge=1 in the the kmer theoretical mass computation
    """
    kmers = generate_kmers(peptide=peptide, max_k=max_kmer_len)
    for kmer in kmers:
        yield Kmer(
            seq=kmer,
            mass=compute_theoretical_mass_over_charge(
                peptide=kmer,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
                charge=charge,
            ),
        )


def get_theoretical_b_y_ion_masses(peptide: str, charge: int = 1):
    b_y_ions = get_b_y_ions(peptide=peptide)
    masses = [
        compute_theoretical_mass_over_charge(peptide=ion, charge=charge)
        for ion in b_y_ions
    ]
    return masses


# def get_b_and_y_ions_of_peptide(peptide: str) -> Tuple[str]:

#     return (compute_theoretical_mass(peptide=aa, charge=charge) for aa in peptide)


# def get_theoretical_b_and_y_ion_mass_spectrum(peptide: str, charge: int = 1) ->


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


# def

# def
