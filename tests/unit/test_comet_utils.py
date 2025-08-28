from src.comet_utils import CometPSM


class Test_CometPSM:
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
    def test_from_assign_confidence_txt(test_data_dir):
        psms = CometPSM.from_txt(
            txt=test_data_dir / "comet_results/assign-confidence.target.txt"
        )
        assert len(psms) > 0
        assert psms[0].q_value is not None
