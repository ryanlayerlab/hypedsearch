import logging
from dataclasses import dataclass
from unittest.mock import patch

import pytest

from src.constants import (
    ION_TYPE_TO_INT,
    MEMORY,
    PRODUCT_ION_TABLE,
    PROTEIN_TABLE,
    IonTypes,
)
from src.peptides_and_ions import Peptide
from src.protein_product_ion_database import (
    DbProductIon,
    DbProtein,
    ProteinProductIonDb,
    create_and_populate_protein_and_product_ion_database,
)
from src.sql_database import PrimaryKey, Sqlite3Database, SqlTableRow
from src.utils import setup_logger
from tests.fixtures_and_helpers import B_NEUTRAL_MASS_CALCULATOR


class Test_DbProtein_from_peptide:
    @staticmethod
    def test_smoke():
        p_id, p_seq = 1, "ABC"
        peptide = Peptide(id=p_id, seq=p_seq)
        actual = DbProtein.from_peptide(peptide=peptide)
        assert actual.id == p_id
        assert actual.seq == p_seq


class Test_DbProductIon_from_peptide:
    @staticmethod
    def test_one_charge_one_ion_type():
        p_seq, p_id = "ACDE", 1
        min_k, max_k = 1, 3
        charge, ion_type = 2, IonTypes.B_ION_TYPE
        peptide = Peptide(seq=p_seq, id=p_id)
        masses = list(range(1, 10))
        with patch(
            B_NEUTRAL_MASS_CALCULATOR,
            side_effect=masses,
        ):
            ions = DbProductIon.from_peptide(
                peptide=peptide,
                charges=[charge],
                ion_types=[ion_type],
                max_k=max_k,
                min_k=min_k,
            )
            assert len(ions) == len(peptide.kmers(min_k=min_k, max_k=max_k))
            assert all([ion.protein_id == p_id for ion in ions])
            assert all([ion.charge == charge for ion in ions])
            assert all(
                [
                    ion.ion_type == ION_TYPE_TO_INT[IonTypes.B_ION_TYPE.value]
                    for ion in ions
                ]
            )
            assert [ion.neutral_mass for ion in ions] == masses

    @staticmethod
    def test_multiple_charges_and_ion_types():
        p_seq, p_id = "ACDE", 1
        min_k, max_k = 1, 3
        charges, ion_types = [1, 2], [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
        peptide = Peptide(seq=p_seq, id=p_id)
        num_kmers = len(peptide.kmers(min_k=min_k, max_k=max_k))
        expected_num_ions = num_kmers * len(charges) * len(ion_types)
        ions = DbProductIon.from_peptide(
            peptide=peptide,
            charges=charges,
            ion_types=ion_types,
            max_k=max_k,
            min_k=min_k,
        )

        assert len(ions) == expected_num_ions
        # Assert number of b-ions and y-ions is half the number of ions
        assert sum(
            [ion.ion_type == ION_TYPE_TO_INT[IonTypes.B_ION_TYPE.value] for ion in ions]
        ) == int(len(ions) / 2)
        assert sum(
            [ion.ion_type == ION_TYPE_TO_INT[IonTypes.Y_ION_TYPE.value] for ion in ions]
        ) == int(len(ions) / 2)

        # Assert number of charge=1 ions and charge=2 ions are both half the number of ions
        assert sum([ion.charge == 1 for ion in ions]) == int(len(ions) / 2)
        assert sum([ion.charge == 2 for ion in ions]) == int(len(ions) / 2)


@dataclass
class PositionTest:
    start: int
    stop: int


@dataclass
class ProteinTest(SqlTableRow):
    id: int
    seq: str


@dataclass
class ProductIonTest(SqlTableRow):
    protein_id: int
    inclusive_start: int
    exclusive_end: int
    charge: int
    neutral_mass: float
    ion_type: int


class Test_ProteinProductIonDb_initialization:
    @staticmethod
    def test_smoke():
        protein_table_name, product_ion_table_name = "proteins", "product_ions"
        min_k, max_k = 2, 4
        charges = [1, 2]
        ion_types = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
        db = ProteinProductIonDb(
            db_path=MEMORY,
            protein_obj=ProteinTest,
            product_ion_obj=ProductIonTest,
            protein_table_name=protein_table_name,
            product_ion_table_name=product_ion_table_name,
            min_k=min_k,
            max_k=max_k,
            charges=charges,
            ion_types=ion_types,
        )

        assert db.db_path == MEMORY
        assert db.charges == charges
        assert db.ion_types == ion_types
        assert db.min_k == min_k
        assert db.max_k == max_k
        assert db.product_ion_obj == ProductIonTest
        assert db.protein_obj == ProteinTest
        assert db.product_ion_table_name == product_ion_table_name
        assert db.protein_table_name == protein_table_name


class Test_ProteinProductionIonDb_add_peptide_and_product_ions:
    @staticmethod
    def test_smoke():

        min_k, max_k = 1, 2
        peptide = Peptide(id=0, seq="ACDE")
        charges = [1]
        ion_types = [IonTypes.B_ION_TYPE]
        db = ProteinProductIonDb(
            min_k=min_k,
            max_k=max_k,
            charges=charges,
            ion_types=ion_types,
        )
        db.create_protein_table()
        db.create_product_ion_table()
        db.add_peptide_and_product_ions(peptide=peptide)
        assert len(db.tables()) == 2

        protein = db.get_proteins()[0]
        assert protein.id == 0
        assert protein.seq == "ACDE"

        product_ions = db.get_product_ions()
        pass
        # assert


class Test_create_and_populate_protein_and_product_ion_database:
    @staticmethod
    def test_smoke(caplog):
        setup_logger()
        p1_seq = "ACDEDE"
        p2_seq = "LMAC"
        peptides = [Peptide(id=0, seq=p1_seq), Peptide(id=1, seq=p2_seq)]
        charges, ion_types = [1], [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
        with caplog.at_level(logging.INFO):
            db = create_and_populate_protein_and_product_ion_database(
                charges=charges,
                ion_types=ion_types,
                peptides=peptides,
            )
        assert len(db.get_proteins()) == 2
        assert len(db.get_product_ions()) > 0
        assert "Adding protein 1 / 2" in caplog.text
        assert "Adding protein 2 / 2" in caplog.text


class Test_ProteinProductIonDb_get_ions_within_mass_tolerance:
    @staticmethod
    def test_smoke():
        masses = [1.8, 1.96, 2.02, 2.5]
        query_mass, mz_tolerance = 2, 0.05
        p_seq = "ACDEF"
        peptides = [Peptide(id=0, seq=p_seq)]
        charges, ion_types = [1], [IonTypes.B_ION_TYPE]

        with patch(
            B_NEUTRAL_MASS_CALCULATOR,
            side_effect=masses,
        ):
            db = create_and_populate_protein_and_product_ion_database(
                charges=charges,
                ion_types=ion_types,
                peptides=peptides,
                min_k=2,
                max_k=2,
            )
            matching_ions = db.get_ions_within_mass_tolerance(
                query_mass=query_mass, mz_tolerance=mz_tolerance
            )

            assert len(matching_ions) == 2
            assert matching_ions[0].neutral_mass == 1.96
            assert matching_ions[1].neutral_mass == 2.02
