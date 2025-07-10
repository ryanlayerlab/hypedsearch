import pytest

from src.constants import GIT_REPO_DIR
from src.hypedsearch_results import (
    HybridPSM,
    NativePSM,
    hypedsearch_output_file_names,
    process_hypedsearch_results,
)
from src.utils import decompress_and_depickle


@pytest.fixture
def hs_result(test_data_dir):
    return decompress_and_depickle(test_data_dir / "hs_SAMPLE_1_878108.pklz")


class Test_process_hypedsearch_results:
    @staticmethod
    def test_pklz_output(tmp_path, hs_result):
        # Arrange
        output_type = "pklz"
        # Act
        process_hypedsearch_results(
            hs_result=hs_result,
            output_type=output_type,
            output_dir=tmp_path,
        )
        output_path = hypedsearch_output_file_names(
            output_dir=tmp_path,
            output_type=output_type,
            mzml_name=hs_result.spectrum.mzml.stem,
            scan=hs_result.spectrum.scan,
        )[0]
        assert output_path.exists()

    @staticmethod
    def test_csv_output(tmp_path, hs_result):
        # Arrange
        output_type = "csv"
        # Act
        process_hypedsearch_results(
            hs_result=hs_result,
            output_type=output_type,
            output_dir=tmp_path,
        )
        output_paths = hypedsearch_output_file_names(
            output_dir=tmp_path,
            output_type=output_type,
            mzml_name=hs_result.spectrum.mzml.stem,
            scan=hs_result.spectrum.scan,
        )
        assert output_paths[0].exists()
        assert output_paths[1].exists()


class Test_NativePSM:
    @staticmethod
    def test_from_comet_psm_and_spectrum(hs_result):
        native_psm = NativePSM.from_native_comet_psm_and_spectrum(
            comet_psm=hs_result.native_psms[0], spectrum=hs_result.spectrum
        )

    @staticmethod
    def test_get_protein_count(hs_result, test_data_dir):
        # Arrange
        protein_counts_path = test_data_dir / "protein_counts.pklz"
        native_psm = NativePSM.from_native_comet_psm_and_spectrum(
            comet_psm=hs_result.native_psms[0], spectrum=hs_result.spectrum
        )
        # Act
        protein_count = native_psm.get_max_protein_count(
            protein_counts=protein_counts_path
        )
        # Assert
        assert protein_count == 2


class Test_HybridPSM:
    @staticmethod
    def test_smoke(hs_result):
        hybrid_psm = HybridPSM.from_hybrid_comet_psm_and_spectrum(
            comet_psm=hs_result.hybrid_psms[0],
            spectrum=hs_result.spectrum,
            hybrid_seqs=hs_result.seq_to_hybrid_peptides,
        )

    @staticmethod
    def test_get_right_side_count(hs_result):
        # Arrange
        comet_results_dir = GIT_REPO_DIR / "results/human_samples/sample1/prot_ab"
        hybrid_psm = HybridPSM.from_hybrid_comet_psm_and_spectrum(
            comet_psm=hs_result.hybrid_psms[0],
            spectrum=hs_result.spectrum,
            hybrid_seqs=hs_result.seq_to_hybrid_peptides,
        )[0]
        # Act
        count = hybrid_psm.get_right_side_ocunt(
            comet_results_dir=comet_results_dir,
        )
        # Assert
        count
