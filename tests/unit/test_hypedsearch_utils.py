import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml

from src.hypedsearch_utils import (
    Crux,
    HypedsearchConfig,
    HypedsearchOutput,
    HypedsearchRunner,
    get_missing_hypedsearch_outputs,
)
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
            runner = Crux()
            assert isinstance(runner, Crux)

        @staticmethod
        def test_crux_exists():
            runner = Crux()

    class Test_run_comet:
        @staticmethod
        def test_successful_run_of_comet(
            tmp_path, mouse_mzml, mouse_fasta, crux_comet_params
        ):
            Crux().run_comet(
                mzml=mouse_mzml,
                fasta=mouse_fasta,
                crux_comet_params=crux_comet_params,
                out_dir=tmp_path,
            )

        @staticmethod
        def test_run_on_spectrum_that_doesnt_produce_psms(
            tmp_path, mouse_mzml, mouse_fasta, crux_comet_params
        ):
            Crux().run_comet(
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
                Crux().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                )
                comet_outputs = Crux.comet_outputs(out_dir=tmp_path)
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=0, fileroot set
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                Crux().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    file_root=file_root,
                )
                comet_outputs = Crux.comet_outputs(
                    out_dir=tmp_path, file_root=file_root
                )
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=1, fileroot unset
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                Crux().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    decoy_search=1,
                )
                comet_outputs = Crux.comet_outputs(out_dir=tmp_path, decoy_search=1)
                assert comet_outputs.target.exists()
                assert comet_outputs.log.exists()
                assert comet_outputs.params.exists()
                assert comet_outputs.decoy is None

            # decoy_search=1, fileroot set
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                Crux().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    file_root=file_root,
                    decoy_search=1,
                )
                comet_outputs = Crux.comet_outputs(
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
                Crux().run_comet(
                    mzml=mouse_mzml,
                    fasta=mouse_fasta,
                    crux_comet_params=crux_comet_params,
                    out_dir=tmp_path,
                    decoy_search=2,
                    scan_min=scan_min,
                    scan_max=scan_max,
                    file_root=file_root,
                )
                comet_outputs = Crux.comet_outputs(
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


class Test_HypedsearchOutput:
    class Test_get_info_from_hypedsearch_comet_output_file:
        @staticmethod
        def test_get_mzml():
            mzml_name = HypedsearchOutput.get_info_from_hypedsearch_comet_output_file(
                file_name="hybrid_mouse_spectra.comet.10-10.decoy.txt", info="mzml"
            )
            assert mzml_name == "mouse_spectra"
            mzml_name = HypedsearchOutput.get_info_from_hypedsearch_comet_output_file(
                file_name="native_mouse_spectra.012.comet.10-10.target.txt", info="mzml"
            )
            assert mzml_name == "mouse_spectra.012"


class Test_HypedsearchConfig:
    class Test_from_yaml:
        @staticmethod
        def test_smoke(test_data_dir):
            # Make sure this passes pydantic's type validation
            HypedsearchConfig.from_yaml(test_data_dir / "hypedsearch_config.yaml")

    class Test_to_yaml:
        @staticmethod
        def test_to_yaml(test_data_dir, tmp_path):
            config = HypedsearchConfig.from_yaml(
                test_data_dir / "hypedsearch_config.yaml"
            )
            out_yaml = tmp_path / "out.yaml"
            config.to_yaml(out_yaml)
            saved_config = HypedsearchConfig.from_yaml(out_yaml)
            assert saved_config == config


class Test_get_missing_hypedsearch_outputs:
    @staticmethod
    def test_smoke():
        config = Path("3A.test.yaml")
        folder = Path("tmp/3A_test")
        missing_outputs = get_missing_hypedsearch_outputs(
            results_dir=folder,
            config=config,
        )
        pass
