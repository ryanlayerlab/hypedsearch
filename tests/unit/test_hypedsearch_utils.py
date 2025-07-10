from pathlib import Path

import pytest

from src.hypedsearch_utils import postprocess_hybrids
from src.mass_spectra import Spectrum
from src.run_hypedsearch import HybridFormingMethod


class Test_postprocess_hybrids:
    @staticmethod
    def test_smoke(test_data_dir):
        mzml = test_data_dir / "BMEM_AspN_Fxn4.mzML"
        db_path = test_data_dir / "top_10_prots.db"
        fasta_path = Path("fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta").absolute()
        peak_to_ion_ppm_tol = 10
        precursor_mz_ppm_tol = 20
        spectra = Spectrum.parse_ms2_from_mzml(mzml)
        hybrid_former = HybridFormingMethod.from_hybrid_formation_relevant_params(
            db_path=db_path,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        )
        hybrids = hybrid_former.form_hybrids(spectrum=spectra[0])
        hybrid_seqs = postprocess_hybrids(hybrids=hybrids, fasta_path=fasta_path)
