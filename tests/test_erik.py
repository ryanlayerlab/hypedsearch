import time
from pathlib import Path
from typing import List

import numpy as np
import pytest
from pyteomics import mzml

from src.erik import (  # Peptide,
    generate_kmers,
    get_b_ion_sequences,
    get_b_y_ion_sequences,
    get_data_for_spectrum,
    get_indices_of_largest_elements,
    get_y_ion_sequences,
    parse_mzml,
    parse_spectrum,
    query_database,
)
from src.erik_constants import (
    CHARGE,
    END,
    ION,
    MASS,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SEQ,
    SPECTRA_DIR,
    START,
    TEST_DIR,
    THOMAS_SAMPLES,
)
from src.fasta_utils import get_proteins_from_fasta
from src.lookups.constants import (
    DOUBLY_CHARGED_B_BASE,
    DOUBLY_CHARGED_Y_BASE,
    SINGLY_CHARGED_B_BASE,
    SINGLY_CHARGED_Y_BASE,
)
from src.lookups.data_classes import Kmer, Peak, Protein, Spectrum
from src.lookups.protein_product_ion_db import ProteinProductIonDb
from tests.fixtures_and_helpers import create_fasta

# from src.runner import (
#     Spectrum,
#     parse_spectrum,
#     relative_abundance_filtering,
#     select_highest_abundance_peaks,
# )

# class TestSelectHighestAbundancePeaks:
#     @staticmethod
#     def test_basic():
#         masses, abundances = (1, 2, 3, 4, 5, 6), (10, 4, 11, 5, 8, 3)
#         num_peaks = 3
#         expected_output = ([1, 3, 5], [10, 11, 8])
#         result = select_highest_abundance_peaks(
#             mass_to_charges=masses, abundances=abundances, num_peaks=num_peaks
#         )
#         assert result == expected_output


# class TestRelativeAbundanceFiltering:
#     @staticmethod
#     def test_basic():
#         masses = [1, 2, 3, 4, 5, 6, 7]
#         abundances = [10, 1, 2, 0, 5, 2, 7]
#         cutoff = 0.3
#         expected_result = ([1, 5, 7], [10, 5, 7])
#         result = relative_abundance_filtering(
#             mass_to_charges=masses, abundances=abundances, rel_abundance_cutoff=cutoff
#         )
#         assert result == expected_result


# class TestParseSpectrum:
#     @staticmethod
#     def test_parse_spectrum():
#         """
#         Test that functions works on all test data we have
#         """
#         # sample = "BMEM_AspN_Fxn5"
#         mzml_paths = list(SPECTRA_DIR.glob("*.mzML"))

#         for file_num, mzml_path in enumerate(mzml_paths):
#             print(f"File {file_num+1}/{len(mzml_paths)}")
#             t0 = time.perf_counter()
#             spectra = [
#                 parse_spectrum(spectrum=spectrum)
#                 for spectrum_num, spectrum in enumerate(mzml.read(str(mzml_path)))
#             ]
#             t1 = time.perf_counter()
#             print(f"Took {t1-t0:.2f} seconds")


def create_spectrum(
    scan_num: int,
    mz_array: List[float] = [1, 2, 3],
    intensity_array: List[float] = [4, 5, 6],
):
    return {
        "id": f"scan={scan_num}",
        "scanList": {"scan": [{"scan start time": 600}]},
        "m/z array": mz_array,
        "intensity array": intensity_array,
        "precursorList": {
            "precursor": [
                {
                    "selectedIonList": {
                        "selectedIon": [
                            {
                                "selected ion m/z": 100,
                                "charge state": 2,
                                "peak intensity": 200,
                            }
                        ]
                    }
                }
            ]
        },
    }


class TestParseSpectrum:
    @staticmethod
    def test_parse_spectrum():
        spectrum = create_spectrum(scan_num=1)
        expected = Spectrum(
            # mass_over_charges=(1, 2, 3),
            # abundances=(4, 5, 6),
            peaks=[
                Peak(mz=1, abundance=4),
                Peak(mz=2, abundance=5),
                Peak(mz=3, abundance=6),
            ],
            precursor_mz=100,
            precursor_charge=2,
            precursor_abundance=200,
            spectrum_id="scan=1",
            retention_time=600,
            scan_num=1,
        )
        actual = parse_spectrum(spectrum)
        assert actual == expected


class TestGetIndicesOfLargestElements:
    @staticmethod
    def test_unique_values():
        array = np.array([5, 8, 1, 4, 6])
        output = get_indices_of_largest_elements(array=array, top_n=3)
        assert (array[output] == [5, 8, 6]).all()

    @staticmethod
    def test_nonunique_values():
        array = np.array([6, 8, 1, 4, 6])
        output = get_indices_of_largest_elements(array=array, top_n=3)
        assert (array[output] == [6, 8, 6]).all()

    @staticmethod
    def test_fewer_values_than_n_value():
        array = np.array([3, 2, 1])
        output = get_indices_of_largest_elements(array=array, top_n=5)
        assert (array[output] == [3, 2, 1]).all()


class TestGenerateKmers:
    @staticmethod
    def test_basic():
        peptide = "ACD"  # B isn't an allowed character
        min_k, max_k = 1, 3
        expected_out = [
            Kmer(seq="A", inclusive_start=0, exclusive_end=1),
            Kmer(seq="C", inclusive_start=1, exclusive_end=2),
            Kmer(seq="D", inclusive_start=2, exclusive_end=3),
            Kmer(seq="AC", inclusive_start=0, exclusive_end=2),
            Kmer(seq="CD", inclusive_start=1, exclusive_end=3),
            Kmer(seq="ACD", inclusive_start=0, exclusive_end=3),
        ]
        actual = generate_kmers(peptide=peptide, min_k=min_k, max_k=max_k)
        assert actual == expected_out


class TestCalculateNeutralMass:
    @staticmethod
    def test_single_charged_y_ion():
        # Arrange
        aa_mass_lookup = {"A": 1, "B": 2, "C": 3}
        seq, charge, ion = "ABBC", 1, "y"
        expected = (
            SINGLY_CHARGED_Y_BASE + sum([aa_mass_lookup[aa] for aa in seq])
        ) / charge

        # Act
        actual = calculate_neutral_mass(
            aa_seq=seq, charge=charge, ion=ion, amino_acid_mass_lookup=aa_mass_lookup
        )

        # Assert
        assert actual == expected

    @staticmethod
    def test_double_charged_y_ion():
        # Arrange
        aa_mass_lookup = {"A": 1, "B": 2, "C": 3}
        seq, charge, ion = "ABBC", 2, "y"
        expected = (
            DOUBLY_CHARGED_Y_BASE + sum([aa_mass_lookup[aa] for aa in seq])
        ) / charge

        # Act
        actual = calculate_neutral_mass(
            aa_seq=seq, charge=charge, ion=ion, amino_acid_mass_lookup=aa_mass_lookup
        )

        # Assert
        assert actual == expected

    @staticmethod
    def test_single_charged_b_ion():
        # Arrange
        aa_mass_lookup = {"A": 1, "B": 2, "C": 3}
        seq, charge, ion = "ABBC", 1, "b"
        expected = (
            SINGLY_CHARGED_B_BASE + sum([aa_mass_lookup[aa] for aa in seq])
        ) / charge

        # Act
        actual = calculate_neutral_mass(
            aa_seq=seq, charge=charge, ion=ion, amino_acid_mass_lookup=aa_mass_lookup
        )

        # Assert
        assert actual == expected

    @staticmethod
    def test_double_charged_b_ion():
        # Arrange
        aa_mass_lookup = {"A": 1, "B": 2, "C": 3}
        seq, charge, ion = "ABBC", 2, "b"
        expected = (
            DOUBLY_CHARGED_B_BASE + sum([aa_mass_lookup[aa] for aa in seq])
        ) / charge

        # Act
        actual = calculate_neutral_mass(
            aa_seq=seq, charge=charge, ion=ion, amino_acid_mass_lookup=aa_mass_lookup
        )

        # Assert
        assert actual == expected


class TestGetBIonSequences:
    @staticmethod
    def test_basic():
        seq = "ABCD"
        expected = ["A", "AB", "ABC", "ABCD"]
        actual = get_b_ion_sequences(peptide="ABCD")
        assert actual == expected


class TestGetYIonSequences:
    @staticmethod
    def test_basic():
        seq = "ABCD"
        expected = ["BCD", "CD", "D"]
        actual = get_y_ion_sequences(peptide="ABCD")
        assert actual == expected


class TestGetBYIonSequences:
    @staticmethod
    def test_basic():
        peptide = "ABCD"
        expected = ["A", "AB", "ABC", "ABCD", "BCD", "CD", "D"]
        actual = get_b_y_ion_sequences(peptide=peptide)
        assert actual == expected


class TestParseMzml:
    @staticmethod
    def test_basic():
        mzml_path = SPECTRA_DIR / f"{THOMAS_SAMPLES[0]}.mzML"
        actual = parse_mzml(mzml_path=mzml_path)
        pass


class TestGetSpecificSpectrum:
    @staticmethod
    def test_basic():
        pass
