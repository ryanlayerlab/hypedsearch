from pathlib import Path

import pytest

from src.mass_spectra import Mzml
from src.process_mzml import (
    aggregate_best_target_psms_and_best_decoy_psms,
    process_mzml,
    process_spectrum,
)


class Test_process_mzml:
    @staticmethod
    def test_via_comet(tmp_path, mouse_fasta, mouse_mzml, mouse_db, comet_params):
        min_scan, max_scan = 1, 30
        process_mzml(
            template_comet_params=comet_params,
            num_psms=20,
            mzml=mouse_mzml,
            fasta=mouse_fasta,
            # out_dir=tmp_path,
            out_dir=Path("tmp/june12/via_comet"),
            db_path=mouse_db,
            min_scan=min_scan,
            max_scan=max_scan,
        )

        assert (tmp_path / "scans_with_psms.json").exists()
        assert (tmp_path / "hybrids_7.json").exists()
        assert (tmp_path / "comet.10-10.target.txt").exists()
        assert (tmp_path / "comet.28-28.decoy.txt").exists()


class Test_process_spectrum:
    @staticmethod
    def test_spectrum_with_hybrids(tmp_path, mouse_fasta, mouse_db):
        mzml = Mzml(path="data/spectra/mouse_samples/BMEM_AspN_Fxn5.mzML")
        scan = 2927
        spectra = [spectrum for spectrum in mzml.ms2_spectra if spectrum.scan == scan]
        assert (
            len(spectra) == 1
        ), f"Expected 1 spectrum for scan {scan}, found {len(spectra)}"
        process_spectrum(
            mzml=mzml.path,
            spectrum=spectra[0],
            db_path=mouse_db,
            precursor_mz_ppm_tol=20.0,
            peak_to_ion_ppm_tol=10.0,
            out_dir=tmp_path,
            num_psms=30,
            fasta=mouse_fasta,
        )

    @staticmethod
    def test_fxn4_scan7(tmp_path, mouse_fasta, mouse_db, comet_params, mouse_mzml):
        mzml = Mzml(path=mouse_mzml)
        scan = 7
        spectra = [spectrum for spectrum in mzml.ms2_spectra if spectrum.scan == scan]
        assert (
            len(spectra) == 1
        ), f"Expected 1 spectrum for scan {scan}, found {len(spectra)}"
        process_spectrum(
            mzml=mzml.path,
            template_comet_params=comet_params,
            spectrum=spectra[0],
            db_path=mouse_db,
            fasta=mouse_fasta,
            out_dir=Path("tmp/june12/via_comet"),
            num_psms=30,
        )


class Test_aggregate_best_target_psms_and_best_decoy_psms:
    @staticmethod
    def test_smoke():
        mzml_dir = Path("tmp/june12/via_comet")
        aggregate_best_target_psms_and_best_decoy_psms(mzml_dir=mzml_dir)

        # assert (mzml_dir / "best_target_psms.json").exists()
        # assert (mzml_dir / "best_decoy_psms.json").exists()
        # assert (mzml_dir / "hybrids.comet.1-10.txt").exists()
