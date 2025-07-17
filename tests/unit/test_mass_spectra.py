from pathlib import Path

import numpy as np

from src.mass_spectra import Mzml, Peak, Spectrum, get_indices_of_largest_elements
from tests.fixtures_and_helpers import create_spectrum


class Test_Mzml:
    @staticmethod
    def test_initialization_by_string(mouse_mzml):
        mzml = Mzml(path=str(mouse_mzml))
        assert isinstance(mzml.path, Path)
        assert mzml.path == mouse_mzml


class Test_Spectrum:
    @staticmethod
    def test_parse_spectrum_from_dict():
        spectrum = create_spectrum(scan_num=1)
        expected = Spectrum(
            peaks=[
                Peak(mz=1, intensity=4, id=0),
                Peak(mz=2, intensity=5, id=1),
                Peak(mz=3, intensity=6, id=2),
            ],
            precursor_mz=100,
            precursor_charge=2,
            precursor_abundance=200,
            spectrum_id="scan=1",
            retention_time=600,
            scan=1,
        )
        actual = Spectrum.from_dict(spectrum=spectrum)
        assert actual == expected

    @staticmethod
    def test_parse_mzml(test_data_dir):
        # It's difficult to create a test MZML. So this test just checks that
        # the function does not fail on an actual MZML
        mzml_path = test_data_dir / "spectra/10_mouse_spectra.mzML"
        actual = Spectrum.parse_ms2_from_mzml(spectra_file=mzml_path)
        assert len(actual) > 0


class Test_get_indices_of_largest_elements:
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
        assert (array[output] == [3, 2, 1]).all()
