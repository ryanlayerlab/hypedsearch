import pandas as pd

from src.constants import B_ION_TYPE, Y_ION_TYPE
from src.hybrids_via_clusters import HybridPeptide
from src.mass_spectra import Spectrum
from src.peptides_and_ions import compute_peptide_precursor_mz
from src.psm import (
    CometPSM,
    PeakIonMatch,
    PeptideSeqSpectrumComparer,
    convert_comet_psms_to_custom_psms,
    get_optimal_ion_support_for_hybrid,
    get_peaks_near_mz,
)


class Test_CometPSM:
    class Test_parsing_txts:
        @staticmethod
        def test_from_comet_txt(comet_txt):
            psms = CometPSM.from_txt(txt=comet_txt)
            assert len(psms) > 0
            assert psms[0].q_value is None

        @staticmethod
        def test_from_crux_txt(crux_txt):
            psms = CometPSM.from_txt(txt=crux_txt)
            assert len(psms) > 0
            assert psms[0].q_value is None

        @staticmethod
        def test_from_assign_confidence_txt(assign_confidence_txt):
            psms = CometPSM.from_txt(txt=assign_confidence_txt)
            assert len(psms) > 0
            assert psms[0].q_value is not None

        @staticmethod
        def test_load_hybrid_psms(test_data_dir):
            psms = CometPSM.from_txt(txt=test_data_dir / "hybrid_comet_psms.txt")
            for psm in psms:
                for prot in psm.proteins:
                    HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)

    class Test_check_if_hybrid_prot:
        @staticmethod
        def test_natives(test_data_dir):
            psms = CometPSM.from_txt(txt=test_data_dir / "native_comet_psms.txt")
            for psm in psms:
                assert not psm.is_hybrid

        @staticmethod
        def test_hybrids(test_data_dir):
            psms = CometPSM.from_txt(txt=test_data_dir / "hybrid_comet_psms.txt")
            for psm in psms:
                assert psm.is_hybrid


class Test_PeptideSeqSpectrumComparer:
    @staticmethod
    def test_init(mouse_mzml):
        comp = PeptideSeqSpectrumComparer(
            spectrum=mouse_mzml.get_spectrum(scan=7),
            seq="SAAPAAGSAPAAAEEKK",
        )
        assert len(list(comp._b_ions_supported)) + len(
            list(comp._y_ions_supported)
        ) == len(comp.peak_ion_matches)
        for b_ion in comp._b_ions_supported:
            assert b_ion.ion_type == B_ION_TYPE
            assert b_ion.ion_name in comp.b_ions_supported_ignore_charge
        for y_ion in comp._y_ions_supported:
            assert y_ion.ion_type == Y_ION_TYPE
            assert y_ion.ion_name in comp.y_ions_supported_ignore_charge

        for num in range(1, 17):
            assert f"b{num}" in comp.ion_name_to_seq_map
            assert len(comp.ion_name_to_seq_map[f"b{num}"]) == num
            assert f"y{num}" in comp.ion_name_to_seq_map
            assert len(comp.ion_name_to_seq_map[f"y{num}"]) == num
        assert comp.num_ions_matched == len(comp.y_ions_supported_with_charge) + len(
            comp.b_ions_supported_with_charge
        )

    @staticmethod
    def test_hybrid_support(mouse_mzml):
        # Arrange
        left_seq = "SAAPAAGSAP"
        right_seq = "AAAEEKK"
        comp = PeptideSeqSpectrumComparer(
            spectrum=mouse_mzml.get_spectrum(scan=7),
            seq="SAAPAAGSAPAAAEEKK",
        )
        # Act
        left_support, right_support = comp.hybrid_support(
            left_seq=left_seq, right_seq=right_seq
        )
        # Assert

    @staticmethod
    def test_peak_to_ion_ppm_differences(mouse_mzml):
        comp = PeptideSeqSpectrumComparer(
            spectrum=mouse_mzml.get_spectrum(scan=7),
            seq="SAAPAAGSAPAAAEEKK",
        )
        comp.peak_ion_matches[0].ppm_diff


class Test_get_optimal_ion_support_for_hybrid:
    @staticmethod
    def test_smoke():
        # Arrange
        expected_left_ions = {"b1", "b2", "b3", "y3", "y4", "y5"}
        expected_right_ions = {"b3", "b4", "b5", "y1", "y2", "y3"}
        # Act
        left_ions, right_ions = get_optimal_ion_support_for_hybrid(
            left_seq="ABC", right_seq="XYZ"
        )
        # Assert
        assert left_ions == expected_left_ions
        assert right_ions == expected_right_ions


class Test_convert_comet_psms_to_custom_psms:
    @staticmethod
    def test_non_hybrids(test_data_dir):
        # Arrange
        psms = CometPSM.from_txt(
            txt=test_data_dir / "BMEM_AspN_Fxn4.assign-confidence.txt"
        )[:10]
        spectra = Spectrum.parse_ms2_from_mzml(
            mzml=test_data_dir / "BMEM_AspN_Fxn4.mzML"
        )
        expected_columns = ["peak_to_ion_mz_ppm_diffs", "precursor_mz_ppm_diff"]
        # Act
        results = convert_comet_psms_to_custom_psms(
            comet_psms=psms,
            spectra=spectra,
        )
        df = pd.DataFrame(results)
        # Assert
        for colm in expected_columns:
            assert colm in df.columns

    @staticmethod
    def test_hybrids(test_data_dir):
        # Arrange
        psms = [
            psm
            for psm in CometPSM.from_txt(
                txt=test_data_dir / "BMEM_AspN_Fxn4.comet.hybrid.target.txt"
            )
            if psm.num == 1
        ][:10]
        spectra = Spectrum.parse_ms2_from_mzml(
            mzml=test_data_dir / "BMEM_AspN_Fxn4.mzML"
        )
        expected_colms = [
            "left_proteins",
            "right_proteins",
            "hybrid_left_support",
            "hybrid_right_support",
            "left_seq",
            "right_seq",
        ]
        # Act
        results = convert_comet_psms_to_custom_psms(
            comet_psms=psms,
            spectra=spectra,
        )
        df = pd.DataFrame(results)
        # Assert
        for colm in expected_colms:
            assert colm in df.columns


class Test_PeakIonMatch:
    @staticmethod
    def test_from_psm(mouse_mzml):
        # Arrange
        spectrum = mouse_mzml.get_spectrum(scan=7)
        seq = "SAAPAAGSAPAAAEEKK"
        ppm_tol = 20
        # Act
        peak_ion_matches = PeakIonMatch.from_psm(
            spectrum=spectrum,
            peptide=seq,
            peak_to_ion_ppm_tolerance=ppm_tol,
        )
        # Assert
        for peak_ion_match in peak_ion_matches:
            assert abs(peak_ion_match.ppm_diff) <= ppm_tol


# def test_smoke(test_data_dir):
#     pass
