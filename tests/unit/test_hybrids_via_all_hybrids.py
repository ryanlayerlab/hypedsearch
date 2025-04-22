import pytest

from src.comet_utils import CometPSM
from src.constants import GIT_REPO_DIR
from src.hybrids_via_all_hybrids import run_hypedsearch_on_spectra

GIT_REPO_DIR


class Test_run_hypedsearch_on_spectrum:
    @staticmethod
    def test_smoke(tmp_path):
        # Arrange
        protein_names = (
            GIT_REPO_DIR / "tests/unit/data/top_10_proteins.txt"
        ).absolute()
        fasta_path = (
            GIT_REPO_DIR / "fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        ).absolute()
        mzml = (GIT_REPO_DIR / "data/spectra/BMEM_AspN_Fxn5.mzML").absolute()
        scan = 9
        output_dir = tmp_path.absolute()
        sample = mzml.stem

        # Act
        run_hypedsearch_on_spectra(
            protein_names=protein_names,
            fasta_path=fasta_path,
            scan=scan,
            mzml=mzml,
            num_psms=20,
            precursor_mz_ppm_tol=20,
            output_dir=output_dir,
        )

        # Assert
        psms = CometPSM.from_txt(file_path=list(output_dir.glob("*"))[0], sample=sample)
        assert len(psms) == 20
