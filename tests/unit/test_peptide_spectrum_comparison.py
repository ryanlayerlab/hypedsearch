import json
from dataclasses import asdict

from src.mass_spectra import Spectrum
from src.peptide_spectrum_comparison import PSM, get_peak_product_ion_matches
from src.utils import mass_difference_in_ppm
from tests.conftest import normalize


class Test_get_peak_product_ion_matches:
    @staticmethod
    def test_smoke(test_data_dir, snapshot, snapshot_dir):
        # Arrange
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        scan = 7
        spectrum = Spectrum.get_spectrum(mzml=mzml, scan=scan)
        peptide = "SAAPAAGSAPAAAEEKK"
        peak_to_ion_ppm_tolerance = 20
        # Act
        peak_ion_matches = get_peak_product_ion_matches(
            spectrum=spectrum,
            peptide=peptide,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
        )
        # Assert
        for peak_ion_match in peak_ion_matches:
            assert (
                mass_difference_in_ppm(
                    mass1=peak_ion_match.peak_mz, mass2=peak_ion_match.ion_mz
                )
                <= peak_to_ion_ppm_tolerance
            )
        peak_ion_matches = normalize([asdict(p) for p in peak_ion_matches])
        peak_ion_matches = sorted(
            peak_ion_matches,
            key=lambda x: (x["ion_seq"], x["ion_charge"], x["peak_mz"]),
        )
        snapshot.snapshot_dir = snapshot_dir
        snapshot_file = f"{spectrum.uid}_peptide={peptide}_peak_to_ion_matches.json"
        snapshot.assert_match(
            json.dumps(peak_ion_matches, indent=2, sort_keys=True),
            snapshot_file,
        )


class Test_PSM:
    @staticmethod
    def test_smoke(test_data_dir):
        spectrum = Spectrum.get_spectrum(
            mzml=test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML", scan=7
        )
        peptide = "SAAPAAGSAPAAAEEKK"
        psm = PSM.from_spectrum_and_seq(
            spectrum=spectrum, seq=peptide, peak_to_ion_ppm_tolerance=20
        )
        # Assert
        assert psm.prefixes_supported == {"SA", "SAA", "SAAPA", "SAAP", "SAAPAA"}
        assert psm.suffixes_supported == {
            "AEEKK",
            "AAPAAGSAPAAAEEKK",
            "APAAGSAPAAAEEKK",
            "AGSAPAAAEEKK",
            "AAGSAPAAAEEKK",
            "GSAPAAAEEKK",
            "APAAAEEKK",
            "AAAEEKK",
            "PAAAEEKK",
            "PAAGSAPAAAEEKK",
            "AAEEKK",
            "SAAPAAGSAPAAAEEKK",
        }
