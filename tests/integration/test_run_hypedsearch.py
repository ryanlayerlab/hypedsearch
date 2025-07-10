from src.mass_spectra import Mzml
from src.run_hypedsearch import process_spectrum


class Test_process_spectrum:
    @staticmethod
    def test_scan_with_native_psms(tmp_path, mouse_mzml, mouse_fasta, mouse_db):
        spectrum = Mzml(path=mouse_mzml).get_spectrum(scan=7)
        process_spectrum(
            spectrum=spectrum,
            db_path=mouse_db,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            num_psms=10,
        )

        expected_out_file_names = [
            "native_target.txt",
            "native_decoy.txt",
            "hybrids.json",
            # "hybrids.fasta",
            "hybrids_target.txt",
            "hybrids_decoy.txt",
        ]
        actual_out_file_names = [out_path.name for out_path in list(tmp_path.glob("*"))]

        for name in expected_out_file_names:
            assert name in actual_out_file_names

    @staticmethod
    def test_scan_with_no_native_psms(tmp_path, mouse_mzml, mouse_fasta, mouse_db):
        spectrum = Mzml(path=mouse_mzml).get_spectrum(scan=1)
        process_spectrum(
            spectrum=spectrum,
            db_path=mouse_db,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            num_psms=10,
        )
        expected_out_file_names = [
            "native_target.txt",
            "native_decoy.txt",
            "hybrids.json",
            # "hybrids.fasta",
            "hybrids_target.txt",
            "hybrids_decoy.txt",
        ]
        actual_out_file_names = [out_path.name for out_path in list(tmp_path.glob("*"))]

        for name in expected_out_file_names:
            assert name in actual_out_file_names
