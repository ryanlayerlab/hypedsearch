from pathlib import Path

import pytest

from src.mass_spectra import (
    get_specific_spectrum_by_sample_and_scan_num,
    get_spectrum_from_mzml,
)
from src.run_hypedsearch import run_hs_on_one_spectrum


class Test_run_hs_on_one_spectrum:
    @staticmethod
    def test_smoke():
        # Arrange
        db_path = Path("results/comet_run_for_protein_abundances/top_10_prots.db")
        mzml_path = Path("data/spectra/BMEM_AspN_Fxn4.mzML")
        scan = 7
        spectrum = get_spectrum_from_mzml(mzml_path=mzml_path, scan_num=scan)

        # Act
        run_hs_on_one_spectrum(db_path=db_path, spectrum=spectrum)
        # Assert

    @staticmethod
    def test_smoke_fxn4_scan7():
        sample = "BMEM_AspN_Fxn4"
        scan_num = 7
        db_path = Path("results/new_fasta/top_10_prots.db").absolute()
        output_dir = Path("tmp")
        spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=sample, scan_num=scan_num
        )
        run_hs_on_one_spectrum(
            db_path=db_path,
            spectrum=spectrum,
            output_dir=output_dir,
            num_peaks=0,
            peak_to_ion_ppm_tol=10,
            precursor_ppm_tol=20,
        )

    @staticmethod
    def test_smoke_fxn5_scan42():
        sample = "BMEM_AspN_Fxn5"
        scan_num = 42
        db_path = Path("results/new_fasta/top_10_prots.db").absolute()
        output_dir = Path("tmp")
        spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=sample, scan_num=scan_num
        )
        run_hs_on_one_spectrum(
            db_path=db_path,
            spectrum=spectrum,
            output_dir=output_dir,
            num_peaks=0,
            peak_to_ion_ppm_tol=10,
            precursor_ppm_tol=20,
        )
