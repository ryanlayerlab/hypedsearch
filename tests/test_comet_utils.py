from pathlib import Path

import pytest

from src.comet_utils import run_comet, update_database_path_in_comet_params_file
from src.constants import COMET_EXECUTABLE, COMET_PARAMS, MOUSE_PROTEOME


class Test_update_database_path_in_comet_params_file:
    @staticmethod
    def test_smoke(tmp_path):
        # Arrange
        fasta_path = "./test.fasta"
        comet_params_path = COMET_PARAMS
        output_path = tmp_path / "test.comet.params"

        # Act
        update_database_path_in_comet_params_file(
            fasta_path=fasta_path,
            comet_params_path=comet_params_path,
            output_path=output_path,
        )

        # Assert
        with open(output_path, "r") as infile:
            new_lines = infile.readlines()
        with open(comet_params_path, "r") as infile:
            old_lines = infile.readlines()
        assert new_lines[4] == "database_name = ./test.fasta\n"
        assert len(new_lines) == len(old_lines)


class Test_run_comet:
    @staticmethod
    def test_keep_params(tmp_path):
        run_comet(
            template_comet_params_path=COMET_PARAMS,
            fasta_path=MOUSE_PROTEOME,
            mzml_path=Path("data/spectra/BMEM_AspN_Fxn4.mzML"),
            comet_exe_path=COMET_EXECUTABLE,
            parent_output_dir=tmp_path,
        )
        comet_params_path = tmp_path / "BMEM_AspN_Fxn4/comet.params"
        comet_results_path = tmp_path / "BMEM_AspN_Fxn4/BMEM_AspN_Fxn4.txt"
        assert comet_params_path.exists()
        assert comet_results_path.exists()
