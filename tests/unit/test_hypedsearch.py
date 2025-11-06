from src.constants import HUMAN_PROTEOME
from src.hypedsearch import TrueHybrid, TrueHybrids
from src.utils import mass_difference_in_ppm


class Test_TrueHybrid:
    @staticmethod
    def test_from_excel():
        path = "data/251022_RP_HUMAN_SERUM_SPIKEDvsNOSPIKE/true_hybrids.xlsx"
        TrueHybrid.from_excel(path=path)

    @staticmethod
    def test_get_proteins_for_true_hybrids():
        path = "data/251022_RP_HUMAN_SERUM_SPIKEDvsNOSPIKE/true_hybrids.xlsx"
        hybrids = TrueHybrid.from_excel(path=path)
        assert hybrids[0].left_proteins == []
        assert hybrids[0].right_proteins == []
        TrueHybrid.get_proteins_for_true_hybrids(
            true_hybrids=hybrids, fasta=HUMAN_PROTEOME
        )
        assert hybrids[0].left_proteins != []
        assert hybrids[0].right_proteins != []

    @staticmethod
    def test_get_spectra_for_true_hybrids():
        # Arrange
        path = "data/251022_RP_HUMAN_SERUM_SPIKEDvsNOSPIKE/true_hybrids.xlsx"
        hybrids = TrueHybrid.from_excel(path=path)
        TrueHybrid.get_proteins_for_true_hybrids(
            true_hybrids=hybrids, fasta=HUMAN_PROTEOME
        )
        mzmls = [
            "data/251022_RP_HUMAN_SERUM_SPIKEDvsNOSPIKE/HumanSerum_with_36_humanHIPS_RP_250430.mzML"
        ]
        precursor_mz_ppm_tol = 20
        rt_tol = 1
        # Act
        TrueHybrid.get_spectra_for_true_hybrids(
            true_hybrids=hybrids,
            mzmls=mzmls,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            retention_time_tol=rt_tol,
        )
        # Assert
        for hybrid in hybrids:
            # Make sure every spectra found is within the chosen precursor m/z tolerance
            # and retention time tolerance
            for spectrum in hybrid.spectra:
                assert (
                    mass_difference_in_ppm(
                        mass1=spectrum.precursor_mz, mass2=hybrid.precursor_mz
                    )
                    <= precursor_mz_ppm_tol
                )
                assert abs(spectrum.retention_time - hybrid.rt)
