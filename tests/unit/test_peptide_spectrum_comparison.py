import json
from dataclasses import asdict

from src.mass_spectra import Spectrum, get_spectrum_from_mzml
from src.peptide_spectrum_comparison import PSM, get_peak_product_ion_matches
from tests.conftest import normalize


class Test_get_peak_product_ion_matches:
    @staticmethod
    def test_smoke(test_data_dir, snapshot, snapshot_dir):
        # Arrange
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        peak_to_ion_ppm_tolerance = 20
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        # Act
        peak_ion_matches = get_peak_product_ion_matches(
            spectrum=spectrum,
            peptide=peptide,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
        )
        # Assert
        peak_ion_matches = normalize([asdict(p) for p in peak_ion_matches])
        peak_ion_matches = sorted(
            peak_ion_matches,
            key=lambda x: (x["ion_seq"], x["ion_charge"], x["peak_mz"]),
        )
        snapshot.snapshot_dir = snapshot_dir
        snapshot_file = f"mzml={mzml.name}_scan={scan_num}_peptide={peptide}_peak_to_ion_matches.json"
        snapshot.assert_match(
            json.dumps(peak_ion_matches, indent=2, sort_keys=True),
            snapshot_file,
        )


class Test_PSM:
    @staticmethod
    def test_init(test_data_dir):
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        psm = PSM(spectrum=spectrum, peptide=peptide)
        assert len(psm.peak_ion_matches) == 20

    @staticmethod
    def test_prefix_support(test_data_dir):
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        psm = PSM(spectrum=spectrum, peptide=peptide)
        assert psm.prefix_support == 5

    @staticmethod
    def test_suffix_support(test_data_dir):
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        psm = PSM(spectrum=spectrum, peptide=peptide)
        assert psm.suffix_support == 12

    @staticmethod
    def test_left_seq_support(test_data_dir):
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        psm = PSM(spectrum=spectrum, peptide=peptide)
        support = psm.left_seq_support("SAAPA")
        assert support == 4

    @staticmethod
    def test_right_seq_support(test_data_dir):
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan_num = 7
        peptide = "SAAPAAGSAPAAAEEKK"
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan_num)
        psm = PSM(spectrum=spectrum, peptide=peptide)
        support = psm.right_seq_support("AEEKK")
        assert support == 1
