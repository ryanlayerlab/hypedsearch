import filecmp
import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
import yaml

from src.constants import COMET_EXECUTABLE, COMET_PARAMS, GIT_REPO_DIR, MOUSE_PROTEOME
from src.mass_spectra import Spectrum
from src.run_comet import (
    CometParamsFile,
    CometRunner,
    get_scans_comet_runs_on,
    run_comet_on_one_spectrum,
    run_comet_via_crux,
)
from src.utils import get_default_comet_executable_path


class Test_CometParams:
    class Test_update_num_output_psms_per_scan:
        @staticmethod
        def test_change_num_psms():
            # Arrange
            num_psms = 20

            # Act
            params = CometParamsFile(num_psms=num_psms)

            # Assert
            for line in params.file_lines:
                if line.strip().startswith("num_output_lines ="):
                    assert line == f"num_output_lines = {num_psms}\n"

    class Test_update_precursor_mz_ppm_tol:
        @staticmethod
        def test_change_precursor_mz_ppm_tol():
            # Arrange
            ppm_tol = 100.0

            # Act
            params = CometParamsFile(precursor_mz_ppm_tol=ppm_tol)

            # Assert
            upper_line_truthiness = False
            lower_line_truthiness = False
            for line in params.file_lines:
                if line.strip().startswith("peptide_mass_tolerance_upper ="):
                    upper_line_truthiness = (
                        line == f"peptide_mass_tolerance_upper = {ppm_tol}\n"
                    )
                if line.strip().startswith("peptide_mass_tolerance_lower ="):
                    lower_line_truthiness = (
                        line == f"peptide_mass_tolerance_lower = -{ppm_tol}\n"
                    )
            assert upper_line_truthiness and lower_line_truthiness

    class Test_update_decoy_search:
        @staticmethod
        def test_update_decoy_search():
            # Arrange
            decoy_search = 2

            # Act
            params = CometParamsFile(decoy_search=decoy_search)

            # Assert
            correct_decoy_search = False
            for line in params.file_lines:
                if line.strip().startswith("decoy_search ="):
                    correct_decoy_search = line == f"decoy_search = {decoy_search}\n"
            assert correct_decoy_search

    class Test_write:
        @staticmethod
        def test_write(tmp_path, snapshot, snapshot_dir):
            # Arrange
            snapshot.snapshot_dir = snapshot_dir
            params = CometParamsFile(
                num_psms=20, precursor_mz_ppm_tol=100.0, decoy_search=1
            )
            output_file = tmp_path / "comet.params"

            # Act
            params.write(output_path=output_file)

            # Assert
            content = output_file.read_text()
            snapshot.assert_match(content, "comet.params")


class Test_CometRunner:
    class Test_comet_output_path:
        @staticmethod
        def test_unset_scan_unset_stem(tmp_path):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            fasta = Path("tests/data/test.fasta")
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=fasta,
                out_dir=tmp_path,
            )

            mzml = Path("blah/test.mzML")
            expected_stem = "test"

            # Act
            out_path = comet_runner._comet_output_path(mzml=mzml)

            # Assert
            assert out_path.stem == expected_stem

        @staticmethod
        def test_set_scan_unset_stem(tmp_path):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            fasta = Path("tests/data/test.fasta")
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=fasta,
                out_dir=tmp_path,
                min_scan=20,
                max_scan=30,
            )

            mzml = Path("blah/test.mzML")
            expected_stem = "test.20-30"

            # Act
            out_path = comet_runner._comet_output_path(mzml=mzml)

            # Assert
            assert out_path.stem == expected_stem

        @staticmethod
        def test_set_scan_set_stem(tmp_path):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            fasta = Path("tests/data/test.fasta")
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=fasta,
                out_dir=tmp_path,
                min_scan=20,
                max_scan=30,
                stem="poop",
            )

            mzml = Path("blah/test.mzML")
            expected_stem = "poop.20-30"

            # Act
            out_path = comet_runner._comet_output_path(mzml=mzml)

            # Assert
            assert out_path.stem == expected_stem

    class Test_run_comet_on_mzml:
        @staticmethod
        def test_output_does_not_exist(tmp_path, mouse_fasta, mouse_mzml):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            fasta = Path("fastas/Uniprot_mouse.fasta")
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=mouse_fasta,
                out_dir=tmp_path,
                max_scan=100,
            )

            # Act
            comet_output_path = comet_runner.run_comet_on_mzml(mzml=mouse_mzml)

            assert comet_output_path.exists()
            assert comet_output_path.stat().st_size > 100

        @staticmethod
        def test_output_exists_overwrite_false(
            caplog, tmp_path, mouse_fasta, mouse_mzml
        ):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=mouse_fasta,
                out_dir=tmp_path,
                stem="poop",
            )
            output_path = tmp_path / "poop.txt"
            output_path.touch()

            # Act
            with caplog.at_level(logging.INFO):
                comet_output_path = comet_runner.run_comet_on_mzml(mzml=mouse_mzml)

            assert "already exist" in caplog.text
            assert comet_output_path.stat().st_size == 0

        @staticmethod
        def test_output_exists_overwrite_true(tmp_path, mouse_fasta, mouse_mzml):
            # Arrange
            comet_exe = get_default_comet_executable_path()
            comet_params = CometParamsFile()
            comet_runner = CometRunner(
                comet_exe=comet_exe,
                comet_params=comet_params,
                fasta=mouse_fasta,
                out_dir=tmp_path,
                stem="poop",
                min_scan=1,
                max_scan=20,
            )
            output_path = tmp_path / "poop.txt"
            output_path.touch()

            # Act
            comet_output_path = comet_runner.run_comet_on_mzml(
                mzml=mouse_mzml, overwrite=True
            )

            assert comet_output_path.stat().st_size > 100


class Test_run_comet_via_crux:
    @staticmethod
    def test_set_scan(tmp_path, mouse_mzml, mouse_fasta):
        out_txt = run_comet_via_crux(
            mzml=mouse_mzml,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            min_scan=1,
            max_scan=10,
            overwrite=True,
        )
        assert out_txt.exists()
        assert out_txt.stat().st_size > 10

    @staticmethod
    def test_all_scans(tmp_path, mouse_mzml, mouse_fasta):
        # Since it's running on all scans in the MZML, this test takes a bit of time
        out_txt = run_comet_via_crux(
            mzml=mouse_mzml,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            overwrite=True,
        )
        assert out_txt.exists()
        assert out_txt.stat().st_size > 10

    @staticmethod
    def test_set_fileroot(tmp_path, mouse_mzml, mouse_fasta):
        # Since it's running on all scans in the MZML, this test takes a bit of time
        out_txt = run_comet_via_crux(
            mzml=mouse_mzml,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            overwrite=True,
            min_scan=1,
            max_scan=10,
            file_root="poop",
        )
        assert out_txt.exists()
        assert out_txt.stem == "poop.comet.1-10"
        assert out_txt.stat().st_size > 10


class Test_get_scans_comet_runs_on:
    @staticmethod
    def test_set_scan_range(tmp_path, mouse_mzml, mouse_fasta):
        scans = get_scans_comet_runs_on(
            out_dir=Path("tmp"),
            fasta=mouse_fasta,
            mzml=mouse_mzml,
            min_scan=1,
            max_scan=10,
        )
        assert len(scans) > 0


class Test_run_comet_on_one_spectrum:
    @staticmethod
    def test_smoke(tmp_path, mouse_mzml):
        comet_params = CometParamsFile()
        spectra = Spectrum.parse_ms2_from_mzml(mouse_mzml)
        scan_num = 7
        spectrum = list(filter(lambda x: x.scan == scan_num, spectra))[0]
        output_path = run_comet_on_one_spectrum(
            template_comet_params=comet_params.template,
            num_psms=comet_params.num_psms,
            precursor_mz_ppm_tolerance=comet_params.precursor_mz_ppm_tol,
            decoy_search=comet_params.decoy_search,
            out_dir=tmp_path,
            fasta=Path("fastas/Uniprot_mouse.fasta"),
            spectrum=spectrum,
            stem="",
            comet_exe=get_default_comet_executable_path(),
            mzml=mouse_mzml,
            overwrite=False,
        )

        assert output_path.exists()
        assert output_path.stat().st_size > 10


class Test_comet_vs_crux_comet:
    @staticmethod
    def test_smoke(mouse_mzml, mouse_fasta, comet_params):
        num_psms = 10
        scan = 7
        crux_output = run_cometcomet_via_crux(
            mzml=mouse_mzml,
            fasta=mouse_fasta,
            out_dir="tmp",
            decoy_search=0,
            num_psms=num_psms,
            min_scan=7,
            max_scan=7,
            file_root="comet_via_crux",
        )

        # run_comet_on_one_spectrum(
        #     template_comet_params=comet_params,
        #     precursor_mz_ppm_tolerance=
        # )
