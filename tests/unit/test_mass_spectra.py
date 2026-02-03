from pathlib import Path

import numpy as np

from src.mass_spectra import Mzml, Peak, Spectrum, get_indices_of_largest_elements
from tests.fixtures_and_helpers import create_spectrum


class Test_Mzml:
    @staticmethod
    def test_initialization_by_string(mouse_mzml_path):
        mzml = Mzml(path=str(mouse_mzml_path))
        assert isinstance(mzml.path, Path)
        assert mzml.path == mouse_mzml_path


class Test_Spectrum:
    @staticmethod
    def test_get_uid():
        sample = "HumanSerum_with_36_humanHIPS_RP_250430"
        scan = 694390
        assert (
            Spectrum.get_uid(sample=sample, scan=scan)
            == "mzml=HumanSerum_with_36_humanHIPS_RP_250430;scan=694390"
        )

    @staticmethod
    def test_parse_uid():
        sample = "HumanSerum_with_36_humanHIPS_RP_250430"
        scan = 694390
        uid = Spectrum.get_uid(sample=sample, scan=scan)
        obs_sample, obs_scan = Spectrum.parse_uid(uid=uid)
        assert sample == obs_sample
        assert scan == obs_scan

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
            precursor_intensity=200,
            spectrum_id="scan=1",
            retention_time=600,
            scan=1,
        )
        actual = Spectrum.from_dict(spectrum=spectrum)
        assert actual == expected

    @staticmethod
    def test_parse_mzml(mouse_mzml_path):
        # It's difficult to create a test MZML. So this test just checks that
        # the function does not fail on an actual MZML
        actual = Spectrum.parse_ms2_from_mzml(mzml=mouse_mzml_path)
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
