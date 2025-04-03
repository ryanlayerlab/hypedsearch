import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from src.comet_utils import CometParams, read_comet_txt_to_df, run_comet_on_one_mzml
from src.constants import COMET_EXECUTABLE, COMET_PARAMS, GIT_REPO_DIR, MOUSE_PROTEOME


class Test_read_comet_txt_to_df:
    @staticmethod
    def test_smoke():
        # Load XML file
        txt_file = GIT_REPO_DIR / "tests/data/test_comet_result.txt"
        df = read_comet_txt_to_df(txt_path=txt_file)
        assert all(df["sample"] == "BMEM_AspN_Fxn5")


class Test_CometParams:
    @staticmethod
    def test_update_database_name(tmp_path):
        # Arrange
        fasta_path = "./test.fasta"
        comet_params_path = COMET_PARAMS

        # Act
        params = CometParams(path=comet_params_path)
        params.update_database_name(fasta_path=fasta_path)

        # Assert
        for line in params.file_lines:
            if line.strip().startswith("database_name ="):
                assert line == f"database_name = {fasta_path}\n"

    @staticmethod
    def test_update_scan_range():
        # Arrange
        min_scan, max_scan = 1, 3
        comet_params_path = COMET_PARAMS

        # Act
        params = CometParams(path=comet_params_path)
        params.update_scan_range(min_scan=min_scan, max_scan=max_scan)

        # Assert
        for line in params.file_lines:
            if line.strip().startswith("scan_range ="):
                assert line == f"scan_range = {min_scan} {max_scan}\n"

    @staticmethod
    def test_num_output_lines():
        # Arrange
        num_output_lines = 100
        comet_params_path = COMET_PARAMS

        # Act
        params = CometParams(path=comet_params_path)
        params.update_num_output_lines(num_output_lines=num_output_lines)

        # Assert
        for line in params.file_lines:
            if line.strip().startswith("num_output_lines ="):
                assert line == f"num_output_lines = {num_output_lines}\n"


class Test_run_comet:
    @staticmethod
    def test_run_with_defaults(tmp_path):
        run_comet_on_one_mzml(
            fasta=MOUSE_PROTEOME,
            mzml=Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute(),
            output_dir=tmp_path,
        )
        comet_params_path = tmp_path / "comet.params"
        comet_results_path = tmp_path / "BMEM_AspN_Fxn4.txt"
        assert comet_params_path.exists()
        assert comet_results_path.exists()

    @staticmethod
    def test_set_stem(tmp_path):
        stem = "poop"
        run_comet_on_one_mzml(
            fasta=MOUSE_PROTEOME,
            mzml=Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute(),
            output_dir=tmp_path,
            stem=stem,
        )
        comet_results_path = tmp_path / f"{stem}.txt"
        assert comet_results_path.exists()

    @staticmethod
    def test_comet_output_exists_and_check_is_false(tmp_path):
        # Arrange
        mzml = Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute()
        comet_output_txt = tmp_path / f"{mzml.stem}.txt"
        comet_output_txt.touch()
        assert comet_output_txt.stat().st_size == 0  # make sure file is empty initially
        # Act
        run_comet_on_one_mzml(
            fasta=MOUSE_PROTEOME,
            mzml=Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute(),
            output_dir=tmp_path,
            check=False,
        )
        # Assert
        assert comet_output_txt.stat().st_size > 0

    @staticmethod
    def test_comet_output_exists_and_check_is_true(tmp_path):
        # Arrange
        mzml = Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute()
        comet_output_txt = tmp_path / f"{mzml.stem}.txt"
        comet_output_txt.touch()
        assert comet_output_txt.stat().st_size == 0  # make sure file is empty initially
        # Act
        run_comet_on_one_mzml(
            fasta=MOUSE_PROTEOME,
            mzml=Path("data/spectra/BMEM_AspN_Fxn4.mzML").absolute(),
            output_dir=tmp_path,
            check=True,
        )
        # Assert
        assert comet_output_txt.stat().st_size == 0
