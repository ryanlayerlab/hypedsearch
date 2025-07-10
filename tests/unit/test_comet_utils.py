import pytest

from src.comet_utils import CometPSM, read_comet_txts_in_dir


class Test_CometPSM:
    @staticmethod
    def test_from_comet_txt(comet_txt):
        psms = CometPSM.from_txt(txt_path=comet_txt)
        assert len(psms) > 0
        assert psms[0].q_value is None

    @staticmethod
    def test_from_crux_txt(crux_txt):
        psms = CometPSM.from_txt(txt_path=crux_txt)
        assert len(psms) > 0
        assert psms[0].q_value is None

    @staticmethod
    def test_from_assign_confidence_txt(test_data_dir):
        psms = CometPSM.from_txt(
            txt_path=test_data_dir / "comet_results/assign-confidence.target.txt"
        )
        assert len(psms) > 0
        assert psms[0].q_value is not None
