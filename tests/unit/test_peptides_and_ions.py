from pathlib import Path
from unittest.mock import patch

from pytest import approx

from src.constants import B_ION_TYPE, Y_ION_TYPE, IonTypes
from src.peptides_and_ions import (
    Fasta,
    Peptide,
    UnpositionedProductIon,
    get_kmer_counts_by_protein,
    get_proteins_by_name,
    get_unique_kmers,
)
from tests.fixtures_and_helpers import (
    B_NEUTRAL_MASS_CALCULATOR,
    Y_NEUTRAL_MASS_CALCULATOR,
)


class Test_ProductIon:
    class Test_compute_b_ion_mz:
        @staticmethod
        def test_default_aa_mass_lookup():
            aa_seq = "ACD"
            charge = 2
            result = UnpositionedProductIon.compute_b_ion_mz(
                seq=aa_seq,
                charge=charge,
            )
            assert result == approx(145.54389711688)

        @staticmethod
        def test_custom_aa_mass_lookup():
            aa_seq = "ACD"
            charge = 2
            aa_mass_lookup = {"A": 1.0, "C": 2.0, "D": 3.0}
            result = UnpositionedProductIon.compute_b_ion_mz(
                seq=aa_seq,
                charge=charge,
                amino_acid_mass_lookup=aa_mass_lookup,
            )
            assert result == approx(4.0072764668800005)

    class Test_compute_y_ion_mz:
        @staticmethod
        def test_default_aa_mass_lookup():
            aa_seq = "ACD"
            charge = 2
            result = UnpositionedProductIon.compute_y_ion_mz(
                seq=aa_seq,
                charge=charge,
            )
            assert result == approx(154.54917946688)

        @staticmethod
        def test_custom_aa_mass_lookup():
            aa_seq = "ACD"
            charge = 2
            aa_mass_lookup = {"A": 1.0, "C": 2.0, "D": 3.0}
            result = UnpositionedProductIon.compute_y_ion_mz(
                seq=aa_seq,
                charge=charge,
                amino_acid_mass_lookup=aa_mass_lookup,
            )
            assert result == approx(13.01255881688)

    class Test_compute_ion_mz:
        @staticmethod
        def test_b_ion():
            seq = "ACD"
            charge = 2
            result = UnpositionedProductIon.compute_ion_mz(
                seq=seq,
                charge=charge,
                ion_type=B_ION_TYPE,
            )
            assert result == approx(145.54389711688)

        @staticmethod
        def test_y_ion():
            seq = "ACD"
            charge = 2
            result = UnpositionedProductIon.compute_ion_mz(
                seq=seq,
                charge=charge,
                ion_type=Y_ION_TYPE,
            )
            assert result == approx(154.54917946688)

    @staticmethod
    def test_get_b_ion_seqs():
        seq = "ACD"
        result = UnpositionedProductIon.get_b_ion_seqs(seq=seq)
        assert result == ["A", "AC", "ACD"]

    @staticmethod
    def test_get_y_ion_seqs():
        seq = "ACD"
        result = UnpositionedProductIon.get_y_ion_seqs(seq=seq)
        assert result == ["ACD", "CD", "D"]

    class Test_generate_product_ions:
        @staticmethod
        def test_default():
            seq = "AC"
            charges = [1, 2]
            ions = UnpositionedProductIon.generate_product_ions(
                seq=seq, charges=charges
            )
            expected_ions = [
                UnpositionedProductIon(seq="A", charge=1, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=1, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="A", charge=2, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=2, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="C", charge=1, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=1, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="C", charge=2, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=2, ion_type=Y_ION_TYPE),
            ]
            assert len(ions) == len(expected_ions)
            for ion in expected_ions:
                assert ion in ions


class Test_get_proteins_by_name:
    @staticmethod
    def test_from_fasta(test_data_dir: Path):
        protein_names = [
            "sp|P01326|INS2_MOUSE",
            "sp|P12968|IAPP_MOUSE",
            "non-existent-protein",
        ]
        prots = get_proteins_by_name(
            protein_names=protein_names,
            fasta_path=test_data_dir
            / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta",
        )
        assert len(prots) == 2
        assert "sp|P01326|INS2_MOUSE" in {p.name for p in prots}
        assert "sp|P12968|IAPP_MOUSE" in {p.name for p in prots}


class Test_Fasta:
    @staticmethod
    def test_write_fasta(tmp_path: Path):
        # Arrange
        peptides = [
            Peptide(seq="ACDEFGHIK", name="peptide1", desc="blah"),
            Peptide(seq="LMNPQRST", name="peptide2"),
        ]
        fasta_path = tmp_path / "test.fasta"
        # Act
        Fasta.write_fasta(peptides=peptides, out_path=fasta_path)
        # Assert
        with open(fasta_path, "r") as f:
            content = f.read()
        expected_content = ">peptide1 blah\nACDEFGHIK\n>peptide2\nLMNPQRST\n"
        assert content == expected_content


def test_get_kmer_ccounts_by_protein(tmp_path: Path):
    # Arrange
    peptides = [
        Peptide(seq="ACACD", name="protein1", desc="blah"),
        Peptide(seq="ACN", name="protein2"),
    ]
    out_path = tmp_path / "test.fasta"
    Fasta.write_fasta(peptides=peptides, out_path=out_path)
    # Act
    kmer_prot_counts = get_kmer_counts_by_protein(fasta=out_path, k=2)
    # Assert
    assert kmer_prot_counts == {
        "AC": {"protein1": 2, "protein2": 1},
        "CA": {"protein1": 1},
        "CD": {"protein1": 1},
        "CN": {"protein2": 1},
    }


class Test_get_unique_kmers:
    @staticmethod
    def test_from_peptides():
        peptides = [
            Peptide(seq="ACACD"),
            Peptide(seq="ACN"),
        ]
        uniq_kmers = get_unique_kmers(peptides=peptides, min_k=1, max_k=2)
        assert uniq_kmers == {"A", "C", "D", "N", "AC", "CA", "CD", "CN"}

    @staticmethod
    def test_from_fasta(tmp_path):
        peptides = [
            Peptide(seq="ACACD", name="protein1"),
            Peptide(seq="ACN", name="protein2"),
        ]
        out_path = tmp_path / "test.fasta"
        Fasta.write_fasta(peptides=peptides, out_path=out_path)
        uniq_kmers = get_unique_kmers(peptides=out_path, min_k=1, max_k=2)
        assert uniq_kmers == {"A", "C", "D", "N", "AC", "CA", "CD", "CN"}
