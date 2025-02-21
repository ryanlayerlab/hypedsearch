from unittest.mock import patch

from src.constants import B_ION_TYPE, Y_ION_TYPE
from src.peptides_and_ions import (
    BIonCreator,
    Fasta,
    Peptide,
    ProductIon,
    YIonCreator,
    generate_product_ions,
    get_proteins_from_fasta,
)
from tests.fixtures_and_helpers import create_fasta

B_NEUTRAL_MASS_CALCULATOR = "src.peptides_and_ions.compute_b_ion_neutral_mass"
Y_NEUTRAL_MASS_CALCULATOR = "src.peptides_and_ions.compute_y_ion_neutral_mass"


class TestGetProteinsFromFasta:
    @staticmethod
    def test_smoke(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Peptide(seq="ATG", desc="seq1", id=0),
            Peptide(seq="CGT", desc="seq2", id=1),
        ]
        actual = get_proteins_from_fasta(fasta_path=fasta_path)
        assert actual == expected_out

    @staticmethod
    def test_via_fasta_class(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Peptide(seq="ATG", desc="seq1", id=0),
            Peptide(seq="CGT", desc="seq2", id=1),
        ]
        fasta = Fasta(path=fasta_path)
        actual = fasta.proteins()
        assert actual == expected_out


class TestProductIonCreators:
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


class TestProductIon:
    @staticmethod
    def test_b_ion_creation():
        seq, charge = "ACD", 2
        ion_type = B_ION_TYPE
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
        ion_type = Y_ION_TYPE
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


class TestGenerateProductIons:
    @staticmethod
    def test_b_ions():
        seq = "AC"
        charges = [1, 2]
        ion_type = B_ION_TYPE
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
                seq=seq, charges=charges, ion_type=ion_type
            )

            assert product_ions == expected

    @staticmethod
    def test_y_ions():
        seq = "AC"
        charges = [1, 2]
        ion_type = Y_ION_TYPE
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
                seq=seq, charges=charges, ion_type=ion_type
            )

            assert product_ions == expected
