from pathlib import Path

from click.testing import CliRunner
from pytest import approx

from src.constants import B_ION_TYPE, MOUSE_PROTEOME, Y_ION_TYPE
from src.peptides_and_ions import (
    Fasta,
    Fasta2MFMIndex,
    Peptide,
    UnpositionedProductIon,
    cli_create_mfm_index_for_fasta,
    get_kmer_counts_by_protein,
    get_proteins_by_name,
    get_unique_kmers,
)
from src.utils import from_pickle


class Test_UnpositionedProductIon:
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

    class Test_generate_product_ions:
        @staticmethod
        def test_default():
            seq = "ACD"
            charges = [1, 2]
            ions = UnpositionedProductIon.generate_product_ions(
                seq=seq, charges=charges
            )
            expected_ions = [
                UnpositionedProductIon(seq="A", charge=1, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=1, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="A", charge=2, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="AC", charge=2, ion_type=B_ION_TYPE),
                UnpositionedProductIon(seq="D", charge=1, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="CD", charge=1, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="D", charge=2, ion_type=Y_ION_TYPE),
                UnpositionedProductIon(seq="CD", charge=2, ion_type=Y_ION_TYPE),
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
            fasta_path=MOUSE_PROTEOME,
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
        Fasta.write_fasta(peptides=peptides, path=fasta_path)
        # Assert
        with open(fasta_path, "r") as f:
            content = f.read()
        expected_content = ">peptide1 blah\nACDEFGHIK\n>peptide2\nLMNPQRST\n"
        assert content == expected_content


def test_get_kmer_counts_by_protein(tmp_path: Path):
    # Arrange
    peptides = [
        Peptide(seq="ACACD", name="protein1", desc="blah"),
        Peptide(seq="ACN", name="protein2"),
    ]
    out_path = tmp_path / "test.fasta"
    Fasta.write_fasta(peptides=peptides, path=out_path)
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
        Fasta.write_fasta(peptides=peptides, path=out_path)
        uniq_kmers = get_unique_kmers(peptides=out_path, min_k=1, max_k=2)
        assert uniq_kmers == {"A", "C", "D", "N", "AC", "CA", "CD", "CN"}


class Test_FastaFMIndex:
    @staticmethod
    def test_from_fasta(tmp_path):
        # Arrange
        peptides = [
            Peptide(seq="ACACD", name="protein1"),
            Peptide(seq="ACN", name="protein2"),
        ]
        path = tmp_path / "prots.fasta"
        Fasta.write_fasta(peptides=peptides, path=path)
        # Act
        converter = Fasta2MFMIndex(fasta=path)
        mfm = converter.create_mfm_index()
        # Assert
        assert mfm.count("AC") == {0: 2, 1: 1}

    @staticmethod
    def test_creation_via_cli(tmp_path):
        # Arrange
        converter = Fasta2MFMIndex(fasta=MOUSE_PROTEOME)
        seq = "KTQILYVMLWLLCVAFTTFLC"
        prot = "sp|Q99LH2|PTSS1_MOUSE"
        # Act
        runner = CliRunner()
        r1 = runner.invoke(
            cli_create_mfm_index_for_fasta,
            ["--fasta", f"{MOUSE_PROTEOME}", "--out_dir", f"{tmp_path}"],
        )
        # Assert
        mfm = from_pickle(path=tmp_path / converter.mfm_name)
        seq_cnt = mfm.count(pattern=seq)
        assert len(seq_cnt.keys()) == 1
        assert converter.idx_to_protein_map[list(seq_cnt.keys())[0]].name == prot
