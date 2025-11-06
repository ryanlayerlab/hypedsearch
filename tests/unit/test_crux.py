from src.comet_utils import CometPSM
from src.constants import MAC_CRUX_EXECUTABLE
from src.crux import Crux, run_comet_on_custom_seqs
from src.mass_spectra import Spectrum


class Test_Crux:
    class Test_run_comet:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir):
            comet_outputs = Crux(path=MAC_CRUX_EXECUTABLE).run_comet(
                mzml=test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML",
                fasta=test_data_dir
                / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta",
                crux_comet_params=test_data_dir / "crux.comet.params",
                decoy_search=2,
                out_dir=tmp_path,
            )
            target_psms = CometPSM.from_txt(txt=comet_outputs.target)
            assert len(target_psms) > 0

    class Test_run_assign_confidence:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir):
            crux = Crux(path=MAC_CRUX_EXECUTABLE)
            out_path = tmp_path / "assign-confidence.txt"
            crux.run_assign_confidence(
                target_txts=[test_data_dir / "BMEM_AspN_Fxn4.target.txt"],
                out_path=out_path,
            )
            psms = CometPSM.from_txt(txt=out_path)
            assert len(psms) > 0


class Test_run_comet_on_custom_seqs:
    @staticmethod
    def test_smoke(test_data_dir):
        seqs = ["EPVDPNRGLRTL", "SAAPAAGSAPAAAEEKK"]
        mzml = test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML"
        spectra = Spectrum.parse_ms2_from_mzml(mzml=mzml)
        spectrum_to_psms = run_comet_on_custom_seqs(
            seqs=seqs,
            spectra=spectra,
            crux_path=MAC_CRUX_EXECUTABLE,
            comet_params=test_data_dir / "crux.comet.params",
        )
        assert isinstance(
            spectrum_to_psms["mzml=BMEM_AspN_Fxn4_scans1-20;scan=1"], CometPSM
        )
        assert spectrum_to_psms["mzml=BMEM_AspN_Fxn4_scans1-20;scan=4"] is None
