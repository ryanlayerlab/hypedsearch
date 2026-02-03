from pathlib import Path

from click.testing import CliRunner

from src.constants import MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.crux import (
    CometConfig,
    CometOutputs,
    Crux,
    cli_run_comet,
    get_expected_comet_outputs_for_mzml_to_scans,
    run_comet_on_custom_seqs,
)
from src.hypedsearch import HypedsearchRunConfig
from src.mass_spectra import Mzml, Spectrum
from src.psm import CometPSM
from src.utils import load_json


class Test_CometOutputs:
    class Test_crux_comet_outputs:
        @staticmethod
        def test_smoke():
            outputs = CometOutputs.crux_comet_outputs(
                out_dir="tmp", decoy_search=0, file_root="name"
            )
            assert str(outputs.target) == "tmp/name.comet.txt"
            assert outputs.decoy is None
            outputs = CometOutputs.crux_comet_outputs(
                out_dir="tmp", decoy_search=1, scan_min=1, scan_max=100
            )
            assert str(outputs.target) == "tmp/comet.1-100.txt"
            assert outputs.decoy is None
            outputs = CometOutputs.crux_comet_outputs(
                out_dir="tmp",
                decoy_search=2,
                file_root="name",
                scan_min=1,
                scan_max=100,
            )
            assert str(outputs.target) == "tmp/name.comet.1-100.target.txt"
            assert str(outputs.decoy) == "tmp/name.comet.1-100.decoy.txt"


class Test_get_expected_comet_outputs_for_mzml_to_scans:
    @staticmethod
    def test_smoke():
        mzml_to_scans = {
            "mzml1": set([1, 2]),
            "mzml2": set([3]),
        }
        outputs = get_expected_comet_outputs_for_mzml_to_scans(
            mzml_to_scans=mzml_to_scans, out_dir="tmp", decoy_search=2, psm_type="both"
        )
        outputs = [str(o) for o in outputs]
        assert len(outputs) == 6
        assert "tmp/comet.1-1.target.txt" in outputs


class Test_Crux:
    class Test_run_comet:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir, mouse_mzml_path):
            comet_outputs = Crux().run_comet(
                mzml=mouse_mzml_path,
                fasta=MOUSE_PROTEOME,
                crux_comet_params=test_data_dir / "crux.comet.params",
                decoy_search=2,
                out_dir=tmp_path,
            )
            target_psms = CometPSM.from_txt(txt=comet_outputs.target)
            assert len(target_psms) > 0

    class Test_run_assign_confidence:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir):
            # Arrange
            crux = Crux()
            out_path = tmp_path / "assign-confidence.target.txt"
            target_txt = test_data_dir / "BMEM_AspN_Fxn4.target.txt"
            # Act
            crux.run_assign_confidence(
                target_txts=[target_txt],
                out_path=out_path,
            )
            # Assert
            psms = CometPSM.from_txt(txt=out_path)
            assert len(psms) > 0


class Test_run_comet_on_custom_seqs:
    @staticmethod
    def test_smoke(test_data_dir, mouse_mzml_path):
        seqs = ["EPVDPNRGLRTL", "SAAPAAGSAPAAAEEKK"]
        spectra = Spectrum.parse_ms2_from_mzml(mzml=mouse_mzml_path)
        spectrum_to_psms = run_comet_on_custom_seqs(
            seqs=seqs,
            spectra=spectra,
            comet_params=test_data_dir / "crux.comet.params",
        )
        assert isinstance(
            spectrum_to_psms["mzml=BMEM_AspN_Fxn4_scans1-20;scan=1"], CometPSM
        )
        assert spectrum_to_psms["mzml=BMEM_AspN_Fxn4_scans1-20;scan=4"] is None
