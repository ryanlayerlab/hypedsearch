import logging
from collections import Counter

import pytest

from src.constants import GIT_REPO_DIR
from src.hybrids_via_clusters import get_hybrids_via_clusters
from src.mass_spectra import Spectrum
from src.peptides_and_ions import KmerToProteinIdMap


class Test_get_hybrids_via_clusters:
    @staticmethod
    def test_smoke(caplog):
        # Arrange
        logger = logging.getLogger("my_module")
        db_path = (GIT_REPO_DIR / "results/new_fasta/top_10_prots.db").absolute()
        mzml = (GIT_REPO_DIR / "data/spectra/BMEM_AspN_Fxn5.mzML").absolute()
        fasta = (
            GIT_REPO_DIR / "fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        ).absolute()
        scan = 9
        # Read MZML to get a spectrum
        spectra = Spectrum.parse_ms2_from_mzml(spectra_file=mzml)
        spectrum = list(filter(lambda spectrum: spectrum.scan == scan, spectra))
        assert len(spectrum) == 1
        spectrum = spectrum[0]

        # Act
        with caplog.at_level(logging.INFO):
            # kmer_to_prot_map = KmerToProteinIdMap.from_fasta(fasta_path=fasta)
            hybrids = get_hybrids_via_clusters(
                db_path=db_path,
                spectrum=spectrum,
                peak_to_ion_ppm_tol=10,
                precursor_mz_ppm_tol=20,
            )

        # Assert
        hybrid_seqs = [hybrid.seq for hybrid in hybrids]
        Counter(hybrid_seqs)
        pass
