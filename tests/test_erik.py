import time
from pathlib import Path
from typing import List

import numpy as np
import pytest
from pyteomics import mzml

from src.erik import (
    Peptide,
    Spectrum,
    generate_kmers,
    get_b_y_ions,
    get_data_for_spectrum,
    get_indices_of_largest_elements,
    parse_spectrum,
)
from src.erik_constants import SPECTRA_DIR, TEST_DIR
from src.erik_utils import Protein, get_proteins_from_fasta
from src.lookups.sqlite_database import Sqllite_Database

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


def test_parse_spectrum():
    spectrum = {
        "id": "scan=1",
        "scanList": {"scan": [{"scan start time": 600}]},
        "m/z array": [1, 2, 3],
        "intensity array": [4, 5, 6],
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
    expected = Spectrum(
        mass_over_charges=(1, 2, 3),
        abundances=(4, 5, 6),
        precursor_mass=100,
        precursor_charge=2,
        precursor_abundance=200,
        spectrum_id="scan=1",
        retention_time=600,
    )
    actual = parse_spectrum(spectrum)
    assert actual == expected


def test_get_data_for_sepctrum():
    sample = "BMEM_AspN_Fxn4"
    scan_num = 2

    get_data_for_spectrum(sample=sample, scan_num=2)


class TestGetBYIons:
    @staticmethod
    def test_basic():
        peptide = "ABCD"
        expected = ["A", "AB", "ABC", "ABCD", "BCD", "CD", "D"]
        actual = get_b_y_ions(peptide=peptide)
        assert actual == expected


class TestPeptide:

    @staticmethod
    def test_b_y_ions():
        peptide = "ABCD"
        expected = ["A", "AB", "ABC", "ABCD", "BCD", "CD", "D"]
        pep = Peptide(peptide=peptide)
        actual = pep.b_y_ions
        assert actual == expected

    # @staticmethod
    # def test_theoretical_b_y_ion_mass_spectrum():


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


def create_fasta(
    folder: Path,
    file_name: str = "test.fasta",
    seqs: List[str] = ["ATGCGTA", "CGTACGT"],
):
    fasta_path = folder / file_name
    lines = []
    for seq_num, seq in enumerate(seqs):
        lines.append(f">seq{seq_num+1}\n")
        lines.append(f"{seq}\n")
    with fasta_path.open("w") as f:
        for line in lines:
            f.write(line)


class TestGetProteinsFromFasta:
    @staticmethod
    def test_basic(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Protein(seq="ATG", desc="seq1"),
            Protein(seq="CGT", desc="seq2"),
        ]

        actual = list(get_proteins_from_fasta(fasta_path=fasta_path))
        assert actual == expected_out


class TestSqlliteDatabase:
    @staticmethod
    def test_insert_proteins(tmp_path):
        # Arrange
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [(0, "seq1", "ATG"), (1, "seq2", "CGT")]
        # Act
        db = Sqllite_Database(path=":memory:", max_len=3)
        proteins = get_proteins_from_fasta(fasta_path=fasta_path)
        db.insert_proteins(proteins=proteins)

        # Assert
        actual = db.cursor.execute("SELECT * FROM proteins").fetchall()
        assert actual == expected_out


class TestGenerateKmers:
    @staticmethod
    def test_basic():
        peptide = "ABCDE"  # B isn't an allowed character
        max_k = 3
        expected_out = [
            # 1-mers
            "A",
            "C",
            "D",
            "E",
            # 2-mers
            "CD",
            "DE",
            # 3-mers
            "CDE",
        ]
        actual = generate_kmers(peptide=peptide, max_k=max_k)
        assert actual == expected_out
