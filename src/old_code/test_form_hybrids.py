import subprocess


from src.form_hybrids import (
    combine_hybrids_from_scans,
    form_hybrids_for_spectrum,
    load_hybrid_peptides_from_json,
)
from src.mass_spectra import Spectrum
from src.utils import load_json


class Test_form_hybrids_for_spectrum:
    @staticmethod
    def test_pass_spectrum_explicitly(mouse_mzml, mouse_fasta, mouse_db):
        seq_to_hy_pep = form_hybrids_for_spectrum(
            spectrum=Spectrum.get_spectrum(scan=7, mzml=mouse_mzml),
            db_path=mouse_db,
            fasta=mouse_fasta,
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
        )
        assert len(seq_to_hy_pep) > 0

    @staticmethod
    def test_pass_spectrum_implicitly(mouse_mzml, mouse_fasta, mouse_db):
        seq_to_hy_pep = form_hybrids_for_spectrum(
            scan=7,
            mzml=mouse_mzml,
            db_path=mouse_db,
            fasta=mouse_fasta,
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
        )
        assert len(seq_to_hy_pep) > 0

    @staticmethod
    def test_save_hybrids_default_out_path(tmp_path, mouse_mzml, mouse_fasta, mouse_db):
        seq_to_hy_pep = form_hybrids_for_spectrum(
            scan=7,
            mzml=mouse_mzml,
            db_path=mouse_db,
            fasta=mouse_fasta,
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
            out_path=tmp_path,
        )
        seq_to_hy_pep = load_json(tmp_path / "7.json")
        assert len(seq_to_hy_pep) > 0

    @staticmethod
    def test_save_hybrids_custom_out_path(tmp_path, mouse_mzml, mouse_fasta, mouse_db):
        out_path = tmp_path / "hybrids.json"
        seq_to_hy_pep = form_hybrids_for_spectrum(
            scan=7,
            mzml=mouse_mzml,
            db_path=mouse_db,
            fasta=mouse_fasta,
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
            out_path=out_path,
        )
        seq_to_hy_pep = load_json(out_path)
        assert len(seq_to_hy_pep) > 0


class Test_cli_form_hybrids:
    @staticmethod
    def test_one_scan(mouse_mzml, mouse_fasta, mouse_db, tmp_path):
        result = subprocess.run(
            [
                "python",
                "-m",
                "src.form_hybrids",
                "form-hybrids",
                "--mzml",
                str(mouse_mzml),
                "--db_path",
                str(mouse_db),
                "--fasta",
                str(mouse_fasta),
                "--out_dir",
                str(tmp_path),
                "--scan",
                "7",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert (tmp_path / "7.json").exists()

    @staticmethod
    def test_many_scans(mouse_mzml, mouse_fasta, mouse_db, tmp_path):
        result = subprocess.run(
            [
                "python",
                "-m",
                "src.form_hybrids",
                "form-hybrids",
                "--mzml",
                str(mouse_mzml),
                "--db_path",
                str(mouse_db),
                "--fasta",
                str(mouse_fasta),
                "--out_dir",
                str(tmp_path),
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        for num in range(1, 11):
            assert (tmp_path / f"{num}.json").exists()


class Test_load_hybrid_peptides_from_json:
    @staticmethod
    def test_smoke(tmp_path, mouse_mzml, mouse_fasta, mouse_db):
        form_hybrids_for_spectrum(
            scan=7,
            mzml=mouse_mzml,
            db_path=mouse_db,
            fasta=mouse_fasta,
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
            out_path=tmp_path,
        )
        load_hybrid_peptides_from_json(tmp_path / "7.json")


class Test_combine_hybrids_from_scans:
    @staticmethod
    def test_dev(tmp_path, mouse_mzml, mouse_db, mouse_fasta):
        result = subprocess.run(
            [
                "python",
                "-m",
                "src.form_hybrids",
                "form-hybrids",
                "--mzml",
                str(mouse_mzml),
                "--db_path",
                str(mouse_db),
                "--fasta",
                str(mouse_fasta),
                "--out_dir",
                str(tmp_path),
            ],
            capture_output=True,
            text=True,
        )
        hybrids = combine_hybrids_from_scans(hybrids_dir=tmp_path, out_path=None)
