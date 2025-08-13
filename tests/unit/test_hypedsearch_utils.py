import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml

from src.hypedsearch_utils import CometRunner, FormHybridsConfig, HypedsearchRunner
from src.mass_spectra import Spectrum


class Test_CometRunner:
    class Test_post_init:
        @staticmethod
        @patch("subprocess.run")
        def test_mock_crux_exists(mock_run):
            # Create a fake result object with returncode=0
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_run.return_value = mock_result
            runner = CometRunner()
            assert isinstance(runner, CometRunner)

        @staticmethod
        def test_crux_exists():
            runner = CometRunner()

    class Test_run_comet:
        @staticmethod
        def test_successful_run_of_comet(
            tmp_path, mouse_mzml, mouse_fasta, crux_comet_params
        ):
            CometRunner().run_comet(
                mzml=mouse_mzml,
                fasta=mouse_fasta,
                crux_comet_params=crux_comet_params,
                out_dir=tmp_path,
            )

        @staticmethod
        def test_run_on_spectrum_that_doesnt_produce_psms(
            tmp_path, mouse_mzml, mouse_fasta, crux_comet_params
        ):
            CometRunner().run_comet(
                mzml=mouse_mzml,
                fasta=mouse_fasta,
                crux_comet_params=crux_comet_params,
                out_dir=tmp_path,
                scan_min=1,
                scan_max=1,
            )

        @staticmethod
        def test_comet_outputs_for_various_settings(
            mouse_mzml, mouse_fasta, crux_comet_params
        ):
            file_root = "testroot_"
            # decoy_search=0, fileroot unset
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                CometRunner().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                )
                comet_outputs = CometRunner.comet_outputs(out_dir=tmp_path)
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=0, fileroot set
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                CometRunner().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    file_root=file_root,
                )
                comet_outputs = CometRunner.comet_outputs(
                    out_dir=tmp_path, file_root=file_root
                )
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=1, fileroot unset
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                CometRunner().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    decoy_search=1,
                )
                comet_outputs = CometRunner.comet_outputs(
                    out_dir=tmp_path, decoy_search=1
                )
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=1, fileroot set
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                CometRunner().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    file_root=file_root,
                    decoy_search=1,
                )
                comet_outputs = CometRunner.comet_outputs(
                    out_dir=tmp_path, decoy_search=1, file_root=file_root
                )
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=2, fileroot set, scan_range set
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                scan_min, scan_max = 1, 8
                CometRunner().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    decoy_search=2,
                    scan_min=scan_min,
                    scan_max=scan_max,
                    file_root=file_root,
                )
                comet_outputs = CometRunner.comet_outputs(
                    out_dir=tmp_path,
                    decoy_search=2,
                    file_root=file_root,
                    scan_min=scan_min,
                    scan_max=scan_max,
                )
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy.exists()


class Test_run_hypedsearch:
    @staticmethod
    def test_non_empty_psms(tmp_path, mouse_mzml, mouse_fasta, crux_comet_params):
        spectrum = Spectrum.get_spectrum(scan=7, mzml=mouse_mzml)
        _ = HypedsearchRunner.run_hypedsearch(
            spectrum=spectrum,
            database=Path("tests/data/mouse_top_10_proteins.db"),
            crux_comet_params=crux_comet_params,
            out_dir=tmp_path,
            fasta=mouse_fasta,
        )
        assert len(list(tmp_path.glob("*"))) == 5
        for file in tmp_path.glob("*"):
            assert file.stat().st_size > 0

    @staticmethod
    def test_empty_psms(tmp_path, mouse_mzml, mouse_fasta, crux_comet_params):
        spectrum = Spectrum.get_spectrum(scan=1, mzml=mouse_mzml)
        _ = HypedsearchRunner.run_hypedsearch(
            spectrum=spectrum,
            database=Path("tests/data/mouse_top_10_proteins.db"),
            crux_comet_params=crux_comet_params,
            out_dir=tmp_path,
            fasta=mouse_fasta,
        )
        assert len(list(tmp_path.glob("*"))) == 5
        for file in tmp_path.glob("*"):
            if ("target" in file.name) or ("decoy" in file.name):
                assert file.stat().st_size == 0
            else:  # hybrids file
                assert file.stat().st_size > 0
