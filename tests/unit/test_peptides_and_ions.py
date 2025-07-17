from unittest.mock import patch


from src.constants import IonTypes
from src.peptides_and_ions import (
    BIonCreator,
    Peptide,
    ProductIon,
    YIonCreator,
    generate_product_ions,
    get_proteins_by_name,
    get_unique_kmers,
    write_fasta,
)
from tests.fixtures_and_helpers import (
    B_NEUTRAL_MASS_CALCULATOR,
    Y_NEUTRAL_MASS_CALCULATOR,
)


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


class Test_get_unique_kmers:
    @staticmethod
    def test_from_list_of_peptides():
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        peptides = [Peptide(seq=p1_seq), Peptide(seq=p2_seq)]
        k = 2
        expected_kmers = {"AC", "CD", "DE", "ED", "LM", "MA"}
        actual_kmers = get_unique_kmers(peptides=peptides, min_k=k, max_k=k)

        assert actual_kmers == expected_kmers


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
    def test_protein_names_as_path(test_data_dir):
        fasta_path = test_data_dir / "test.fasta"
        protein_names = ["sp|B0V2N1|PTPRS_MOUSE"]
        prots = get_proteins_by_name(protein_names=protein_names, fasta_path=fasta_path)
        assert len(prots) == 1
        assert prots[0].name == protein_names[0]
