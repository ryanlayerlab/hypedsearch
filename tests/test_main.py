from main import hypedsearch
from src.erik_constants import MOUSE_PROTEOME, SPECTRA_DIR, TEST_DIR, THOMAS_SAMPLES


class TestHypedsearch:
    @staticmethod
    def test_full_workflow(tmp_path):
        # Arrange
        mzml_path = SPECTRA_DIR / f"{THOMAS_SAMPLES[0]}.mzML"
        ppm_tol = 5
        fasta_path = TEST_DIR / "test.fasta"
        max_kmer_len = 30
        charges_to_consider = [1]
        db_path = tmp_path / "protein_kmer.db"

        # Act
        hypedsearch(
            ppm_tol=ppm_tol,
            mzml_path=mzml_path,
            fasta_path=fasta_path,
            max_kmer_len=max_kmer_len,
            db_path=db_path,
            charges_to_consider=charges_to_consider,
        )
