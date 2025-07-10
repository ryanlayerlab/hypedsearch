from pathlib import Path
from unittest.mock import MagicMock

import pytest

from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import get_specific_spectrum_by_sample_and_scan_num
from src.native_and_hybrid_comet_run import CometRun, run_hypedsearch

COMET_TXT = Path("comet/comet.macos.exe").absolute()
COMET_PARAMS = Path("comet/comet.params").absolute()
FASTA = Path("fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta").absolute()
assert COMET_TXT.exists()
assert COMET_PARAMS.exists()
assert FASTA.exists()


class Test_CometRun:
    @staticmethod
    def test_initialization():
        CometRun(
            exe=COMET_TXT,
            params=COMET_PARAMS,
            fasta=FASTA,
            num_psms=3,
        )

    @staticmethod
    def test_run():
        # Arrange
        comet_run = CometRun(
            exe=COMET_TXT, params=COMET_PARAMS, fasta=FASTA, num_psms=3
        )
        sample = "BMEM_AspN_Fxn4"
        scan = 7
        spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=sample, scan_num=scan
        )
        # Act
        psms = comet_run.run_on_spectrum(spectrum=spectrum)
        # Assert
        assert len(psms) == 3


class Test_run_hypedsearch:
    @staticmethod
    def test_smoke():
        # Arrange
        sample = "BMEM_AspN_Fxn4"
        scan = 7
        spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=sample, scan_num=scan
        )
        mock_hybrid_former = MagicMock()
        hy_pep = HybridPeptide(b_seq="AAA", y_seq="BBB")
        hy_pep.fasta_description = "b-prots:AAA y-prots:BBB"
        mock_hybrid_former.form_hybrids.return_value = [hy_pep]
        # Act
        hs_result = run_hypedsearch(
            spectrum=spectrum,
            hybrid_former=mock_hybrid_former,
            fasta=FASTA,
            comet_exe=COMET_TXT,
            comet_params=COMET_PARAMS,
            other_prots_for_fasta=FASTA,
            num_psms=3,
        )
        # Assert
        # The hybrid isn't a could match for the spectrum so the PSMs between
        # runs shouldn't change
        assert hs_result.hybrid_psms == hs_result.native_psms
