from src.peptide_spectrum_comparison import PeptideSpectrumComparison
from src.utils import decompress_and_depickle


class Test_PeptideSpectrumComparison:
    @staticmethod
    def test_native_psm():
        # Arrange
        hs_result = decompress_and_depickle(
            "results/human_samples/sample1/hs_results/hs_SAMPLE_1_696578.pklz"
        )
        # Act
        comp = PeptideSpectrumComparison.from_str(
            spectrum=hs_result.spectrum,
            peptide=hs_result.native_psms[0].seq,
        )
        # Assert
        assert comp.num_hybrid_seqs is None
        assert comp.hybrid_seq is None
        assert comp.num_b_ion_seqs == len(hs_result.native_psms[0].seq)
        assert comp.num_y_ion_seqs == len(hs_result.native_psms[0].seq)

    @staticmethod
    def test_hybrid_psm():
        # Arrange
        hs_result = decompress_and_depickle(
            "results/human_samples/sample1/hs_results/hs_SAMPLE_1_696578.pklz"
        )
        hybrid_psm = hs_result.hybrid_psms[0]
        hybrid_peptide = hs_result.hybrid_seqs[hybrid_psm.seq][0]
        # Act
        comp = PeptideSpectrumComparison.from_str(
            spectrum=hs_result.spectrum,
            peptide=hybrid_peptide.seq_with_hyphen,
        )
        # Assert
        # tmp = comp.product_ion_seqs_with_matching_peaks
        assert comp.hybrid_seq == "PQGAWKEL-GVQAPE"
        assert comp.num_hybrid_seqs == len("PQGAWKELGVQAPE")
        assert comp.num_hybrid_seqs_supported == 1
