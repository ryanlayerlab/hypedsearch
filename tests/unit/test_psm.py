from src.hybrids_via_clusters import HybridPeptide
from src.psm import CometPSM


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
