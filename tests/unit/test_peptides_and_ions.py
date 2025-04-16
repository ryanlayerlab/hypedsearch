from unittest.mock import patch

import pytest

from src.constants import B_ION_TYPE, GIT_REPO_DIR, MOUSE_PROTEOME, Y_ION_TYPE, IonTypes
from src.peptides_and_ions import (
    BIonCreator,
    Peptide,
    ProductIon,
    YIonCreator,
    generate_product_ions,
    get_proteins_by_name,
    get_proteins_from_fasta,
    get_uniq_kmer_to_protein_map,
    get_unique_kmers,
    random_sample_of_unique_kmers,
    write_fasta,
)
from tests.fixtures_and_helpers import (
    B_NEUTRAL_MASS_CALCULATOR,
    Y_NEUTRAL_MASS_CALCULATOR,
    create_fasta,
)


class Test_get_proteins_from_fasta:
    @staticmethod
    def test_smoke(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Peptide(seq="ATG", name="seq1", desc="seq1", id=0),
            Peptide(seq="CGT", name="seq2", desc="seq2", id=1),
        ]
        actual = get_proteins_from_fasta(fasta_path=fasta_path)
        assert actual == expected_out

    @staticmethod
    def test_via_fasta_class(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Peptide(seq="ATG", name="seq1", desc="seq1", id=0),
            Peptide(seq="CGT", name="seq2", desc="seq2", id=1),
        ]
        actual = Peptide.from_fasta(fasta_path=fasta_path)
        assert actual == expected_out


class Test_ProductIonCreators_generate_product_ion_seqs:
    @staticmethod
    def test_b_ions():
        seq = "ABCD"
        y_seqs = BIonCreator().generate_product_ion_seqs(seq=seq)
        expected = ["A", "AB", "ABC", "ABCD"]
        assert y_seqs == expected

    @staticmethod
    def test_y_ions():
        seq = "ABCD"
        b_seqs = YIonCreator().generate_product_ion_seqs(seq=seq)
        expected = ["ABCD", "BCD", "CD", "D"]
        assert b_seqs == expected


class Test_ProductIon_initialization:
    @staticmethod
    def test_b_ion_creation():
        seq, charge = "ACD", 2
        ion_type = IonTypes.B_ION_TYPE
        neutral_mass = 1.5

        # Patch the b-ion neutral mass calculation to separate
        # tests for that function from tests of ProductIon instantiation
        with patch(
            B_NEUTRAL_MASS_CALCULATOR,
            return_value=neutral_mass,
        ):
            ion = ProductIon(seq=seq, charge=charge, ion_type=ion_type)

            assert ion.seq == seq
            assert ion.charge == charge
            assert ion.ion_type == ion_type
            assert ion.neutral_mass == neutral_mass

    @staticmethod
    def test_y_ion_creation():
        seq, charge = "ACD", 2
        ion_type = IonTypes.Y_ION_TYPE
        neutral_mass = 2.0

        # Patch the b-ion neutral mass calculation to separate
        # tests for that function from tests of ProductIon instantiation
        with patch(
            Y_NEUTRAL_MASS_CALCULATOR,
            return_value=neutral_mass,
        ):
            ion = ProductIon(seq=seq, charge=charge, ion_type=ion_type)

            assert ion.seq == seq
            assert ion.charge == charge
            assert ion.ion_type == ion_type
            assert ion.neutral_mass == neutral_mass


class Test_generate_product_ions:
    @staticmethod
    def test_b_ions():
        seq = "AC"
        charges = [1, 2]
        ion_type = IonTypes.B_ION_TYPE
        ion_types = [ion_type]
        # Patch neutral mass calculator to decouple tests for this function
        # from tests of neutral mass calculator
        with patch(
            B_NEUTRAL_MASS_CALCULATOR,
            side_effect=[1, 2, 1, 2, 1, 2, 1, 2],
        ):
            expected = [
                ProductIon(
                    seq="A",
                    charge=1,
                    ion_type=ion_type,
                    # neutral_mass will be 1
                ),
                ProductIon(
                    seq="AC",
                    charge=1,
                    ion_type=ion_type,
                    # neutral_mass will be 2
                ),
                ProductIon(
                    seq="A",
                    charge=2,
                    ion_type=ion_type,
                    # neutral_mass will be 1
                ),
                ProductIon(
                    seq="AC",
                    charge=2,
                    ion_type=ion_type,
                    # neutral_mass will be 2
                ),
            ]

            product_ions = generate_product_ions(
                seq=seq, charges=charges, ion_types=ion_types
            )

            assert product_ions == expected

    @staticmethod
    def test_y_ions():
        seq = "AC"
        charges = [1, 2]
        ion_type = IonTypes.Y_ION_TYPE
        ion_types = [ion_type]
        # Patch neutral mass calculator to decouple tests for this function
        # from tests of neutral mass calculator
        with patch(
            Y_NEUTRAL_MASS_CALCULATOR,
            side_effect=[1, 2, 1, 2, 1, 2, 1, 2],
        ):
            expected = [
                ProductIon(
                    seq="AC",
                    charge=1,
                    ion_type=ion_type,
                    # neutral_mass will be 1
                ),
                ProductIon(
                    seq="C",
                    charge=1,
                    ion_type=ion_type,
                    # neutral_mass will be 2
                ),
                ProductIon(
                    seq="AC",
                    charge=2,
                    ion_type=ion_type,
                    # neutral_mass will be 1
                ),
                ProductIon(
                    seq="C",
                    charge=2,
                    ion_type=ion_type,
                    # neutral_mass will be 2
                ),
            ]

            product_ions = generate_product_ions(
                seq=seq, charges=charges, ion_types=ion_types
            )

            assert product_ions == expected


class TestRandomSampleOfUniqueKmers:
    @staticmethod
    def test_one_peptide():
        p_seq = "ACDE"
        peptides = [Peptide(seq=p_seq)]
        k, sample_size = 2, 3
        kmers = random_sample_of_unique_kmers(
            peptides=peptides, k=k, sample_size=sample_size
        )
        assert len(kmers) == sample_size
        assert all([kmer in p_seq for kmer in kmers])

    @staticmethod
    def test_two_peptide():
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        peptides = [Peptide(seq=p1_seq), Peptide(seq=p2_seq)]
        k, sample_size = 2, 5
        kmers = random_sample_of_unique_kmers(
            peptides=peptides, k=k, sample_size=sample_size
        )
        assert len(kmers) == sample_size
        assert all([((kmer in p1_seq) or (kmer in p2_seq)) for kmer in kmers])

    @staticmethod
    def test_raise_error():
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        peptides = [Peptide(seq=p1_seq), Peptide(seq=p2_seq)]
        k, sample_size = 2, 20
        with pytest.raises(RuntimeError, match="number of unique kmers"):
            random_sample_of_unique_kmers(
                peptides=peptides, k=k, sample_size=sample_size
            )

    @staticmethod
    def test_via_fasta(tmp_path):
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        file_name = "test.fasta"
        create_fasta(seqs=[p1_seq, p2_seq], folder=tmp_path, file_name=file_name)
        k, sample_size = 2, 4
        kmers = random_sample_of_unique_kmers(
            fasta_path=tmp_path / file_name, k=k, sample_size=sample_size
        )
        assert len(kmers) == sample_size
        assert all([((kmer in p1_seq) or (kmer in p2_seq)) for kmer in kmers])


class Test_get_unique_kmers:
    @staticmethod
    def test_smoke():
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        peptides = [Peptide(seq=p1_seq), Peptide(seq=p2_seq)]
        k = 2
        expected_kmers = {"AC", "CD", "DE", "ED", "LM", "MA"}
        actual_kmers = get_unique_kmers(peptides=peptides, k=k)

        assert actual_kmers == expected_kmers


class Test_get_uniq_kmer_to_protein_map:
    @staticmethod
    def test_smoke():
        proteins = [Peptide(seq="ACD", id=0), Peptide(seq="CDE", id=1)]
        min_k, max_k = 1, 3
        expected = {
            "A": [0],
            "AC": [0],
            "ACD": [0],
            "C": [0, 1],
            "CD": [0, 1],
            "CDE": [1],
            "D": [0, 1],
            "DE": [1],
            "E": [1],
        }
        actual = get_uniq_kmer_to_protein_map(
            min_k=min_k, max_k=max_k, proteins=proteins
        )

        assert actual == expected


class Test_write_fasta:
    @staticmethod
    def test_smoke(tmp_path):
        output_path = tmp_path / "test.fasta"
        peptides = [
            Peptide("ABC", name="blah1", desc="this is blah1", id=0),
            Peptide("DEF", name="blah2", desc="this is blah2", id=1),
        ]
        expected_lines = [">blah1 this is blah1", "ABC", ">blah2 this is blah2", "DEF"]
        write_fasta(peptides=peptides, output_path=output_path)
        lines = []
        with open(output_path, "r") as fasta:
            lines = [line.strip() for line in fasta]

        assert lines == expected_lines


class Test_get_proteins_by_name:
    @staticmethod
    def test_protein_names_as_path():
        fasta_path = GIT_REPO_DIR / "fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        protein_names = GIT_REPO_DIR / "tests/unit/data/top_10_proteins.txt"
        prots = get_proteins_by_name(protein_names=protein_names, fasta_path=fasta_path)
        assert len(prots) == 10

    # @staticmethod
