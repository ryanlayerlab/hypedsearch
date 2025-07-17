import logging
from dataclasses import dataclass
from typing import Union

import pytest

from src.create_db import KmerDatabase
from src.hybrids_via_clusters import SpectrumExtendedClusters
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Peptide


@dataclass
class PositionedIon:
    seq: str
    charge: int
    ion_type: str
    protein: Union[int, str]
    inclusive_start: int
    exclusive_end: int

    # @classmethod
    # def from_ion_and_location(
    #     cls,
    #     ion: ProductIon,
    # ) -> "PositionedIon":
    #     return cls(
    #         seq=seq,
    #         charge=charge,
    #         ion_type=ion_type.value,
    #         protein=protein,
    #         inclusive_start=inclusive_start,
    #         exclusive_end=exclusive_end,
    #     )


class Test_SpectrumExtendedClusters:
    class Test_from_spectrum:
        @staticmethod
        def test_smoke(test_data_dir):
            scan = 7
            mzml = "data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML"
            proteins = test_data_dir / "test.fasta"
            db_path = test_data_dir / "test.fasta.db"
            peak_to_ion_ppm_tol, precursor_mz_ppm_tol = 20, 20
            spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
            spectrum_extended_clusters = SpectrumExtendedClusters.from_spectrum(
                spectrum=spectrum,
                db_path=db_path,
                proteins=proteins,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            )
            hybrids = spectrum_extended_clusters.form_hybrids(
                precursor_charge=spectrum.precursor_charge,
                precursor_mz=spectrum.precursor_mz,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            )
            assert len(hybrids) > 0

    # @staticmethod
    # def test_testing(test_data_dir):
    #     scan = 7
    #     mzml = "data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML"
    #     fasta = test_data_dir / "test.fasta"
    #     db_path = test_data_dir / "test.fasta.db"
    #     spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)

    #     db = KmerDatabase(db_path=db_path)
    #     # Get peak-ion matches for the spectrum
    #     peak_ion_matches = db.get_spectrum_peak_ion_matches(
    #         spectrum=spectrum,
    #         ppm_tolerance=20,
    #     )
    #     # For each peak-ion match, find all its potential locations in the proteins
    #     protein_name_to_seq_map = {
    #         pep.name: pep.seq for pep in Peptide.from_fasta(fasta)
    #     }
    #     # For each peak-ion match, find all its potential locations in the proteins
    #     positioned_ions = []
    #     for peak_ion_match in peak_ion_matches:
    #         for protein in peak_ion_match.ion.proteins:
    #             locations_in_protein = get_positions_of_subseq_in_seq(
    #                 subseq=peak_ion_match.seq,
    #                 seq=protein_name_to_seq_map[protein],
    #             )
    #             positioned_ions.extend(
    #                 [
    #                     PositionedIon(
    #                         seq=peak_ion_match.seq,
    #                         charge=peak_ion_match.ion.charge,
    #                         ion_type=peak_ion_match.ion.ion_type,
    #                         protein=protein,
    #                         inclusive_start=loc.inclusive_start,
    #                         exclusive_end=loc.exclusive_end,
    #                     )
    #                     for loc in locations_in_protein
    #                 ]
    #             )
    #     pass
