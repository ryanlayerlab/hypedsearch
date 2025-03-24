import pytest

from src.constants import IonTypes
from src.mass_spectra import Peak
from src.peptide_spectrum_comparison import (
    Cluster,
    get_b_clusters,
    get_ions_matching_peak,
    get_y_clusters,
)
from src.protein_product_ion_database import (
    B_ION_TYPE,
    DbKmer,
    DbProtein,
    ProteinProductIonDb,
)


class Test_get_ions_matching_peak:
    @staticmethod
    def test_smoke():
        # Arrange
        min_k, max_k = 1, 3
        charges = [1]
        ion_types = [IonTypes.B_ION_TYPE]
        proteins = [DbProtein(id=1, seq="ABC")]
        product_ions = [
            DbKmer(
                protein_id=1,
                inclusive_start=1,
                exclusive_end=2,
                charge=1,
                neutral_mass=1.2,
                ion_type=0,
            ),
            DbKmer(
                protein_id=1,
                inclusive_start=1,
                exclusive_end=3,
                charge=1,
                neutral_mass=0.9,
                ion_type=0,
            ),
            DbKmer(
                protein_id=1,
                inclusive_start=2,
                exclusive_end=3,
                charge=1,
                neutral_mass=3,
                ion_type=0,
            ),
            # 2-mers
            DbKmer(
                protein_id=1,
                inclusive_start=0,
                exclusive_end=2,
                charge=1,
                neutral_mass=4,
                ion_type=0,
            ),
        ]
        peak = Peak(mz=1, intensity=100)
        # Act
        db = ProteinProductIonDb(
            min_k=min_k,
            max_k=max_k,
            charges=charges,
            ion_types=ion_types,
        )
        db.create_protein_table()
        db.create_product_ion_table()
        db.insert_dataclasses(table_name=db.protein_table_name, data_classes=proteins)
        db.insert_dataclasses(
            table_name=db.product_ion_table_name, data_classes=product_ions
        )

        peak_with_matches = get_ions_matching_peak(
            peak=peak, precursor_charge=2, ppm_tolerance=300000, db=db
        )
        assert len(peak_with_matches.ions) == 2
        assert peak_with_matches.ions == product_ions[0:2]


class Test_Cluster_initialization:
    @staticmethod
    def test_b_cluster():
        cluster = Cluster(
            ions=[
                DbKmer(
                    protein_id=0,
                    inclusive_start=1,
                    exclusive_end=3,
                    charge=2,
                    neutral_mass=5.2,
                    ion_type=0,
                ),
                DbKmer(
                    protein_id=0,
                    inclusive_start=1,
                    exclusive_end=5,
                    charge=2,
                    neutral_mass=5.2,
                    ion_type=0,
                ),
            ]
        )
        assert cluster.protein_id == 0
        assert cluster.inclusive_start == 1
        assert cluster.exclusive_end == 5
        assert cluster.ion_type == "b"

    @staticmethod
    def test_y_cluster():
        cluster = Cluster(
            ions=[
                DbKmer(
                    protein_id=1,
                    inclusive_start=1,
                    exclusive_end=5,
                    charge=2,
                    neutral_mass=5.2,
                    ion_type=1,
                ),
                DbKmer(
                    protein_id=1,
                    inclusive_start=3,
                    exclusive_end=5,
                    charge=2,
                    neutral_mass=5.2,
                    ion_type=1,
                ),
            ]
        )
        assert cluster.protein_id == 1
        assert cluster.inclusive_start == 1
        assert cluster.exclusive_end == 5
        assert cluster.ion_type == "y"

    @staticmethod
    def test_inconsistient_ion_type():
        with pytest.raises(
            AssertionError, match="All ions should have the same ion type!"
        ):
            Cluster(
                ions=[
                    DbKmer(
                        protein_id=0,
                        inclusive_start=1,
                        exclusive_end=3,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=0,
                    ),
                    DbKmer(
                        protein_id=0,
                        inclusive_start=1,
                        exclusive_end=5,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=1,
                    ),
                ]
            )

    @staticmethod
    def test_inconsistient_start_locations():
        with pytest.raises(
            AssertionError, match="All b-ions should have same start position"
        ):
            Cluster(
                ions=[
                    DbKmer(
                        protein_id=0,
                        inclusive_start=1,
                        exclusive_end=3,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=0,
                    ),
                    DbKmer(
                        protein_id=0,
                        inclusive_start=2,
                        exclusive_end=5,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=0,
                    ),
                ]
            )

    @staticmethod
    def test_inconsistient_end_locations():
        with pytest.raises(
            AssertionError, match="All y-ions should have same end position"
        ):
            Cluster(
                ions=[
                    DbKmer(
                        protein_id=0,
                        inclusive_start=1,
                        exclusive_end=3,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=1,
                    ),
                    DbKmer(
                        protein_id=0,
                        inclusive_start=2,
                        exclusive_end=5,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=1,
                    ),
                ]
            )

    @staticmethod
    def test_inconsistient_protein_ids():
        with pytest.raises(
            AssertionError, match="All ions should have same protein ID"
        ):
            Cluster(
                ions=[
                    DbKmer(
                        protein_id=0,
                        inclusive_start=1,
                        exclusive_end=3,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=0,
                    ),
                    DbKmer(
                        protein_id=1,
                        inclusive_start=1,
                        exclusive_end=5,
                        charge=2,
                        neutral_mass=5.2,
                        ion_type=0,
                    ),
                ]
            )


class Test_get_b_clusters:
    @staticmethod
    def test_smoke():
        # Arrange
        ions = [
            # b-cluster 1
            DbKmer(
                protein_id=0,
                inclusive_start=1,
                exclusive_end=3,
                charge=2,
                neutral_mass=5.2,
                ion_type=0,
            ),
            DbKmer(
                protein_id=0,
                inclusive_start=1,
                exclusive_end=5,
                charge=2,
                neutral_mass=5.2,
                ion_type=0,
            ),
            # b-cluster 2
            DbKmer(
                protein_id=1,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=0,
            ),
            # y-ions
            DbKmer(
                protein_id=1,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=1,
            ),
            DbKmer(
                protein_id=3,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=1,
            ),
        ]

        # Act
        b_clusters = get_b_clusters(ions=ions)

        # Assert
        assert len(b_clusters) == 2

        assert b_clusters[0].protein_id == 0
        assert len(b_clusters[0].ions) == 2
        assert b_clusters[0].inclusive_start == 1
        assert b_clusters[0].exclusive_end == 5

        assert b_clusters[1].protein_id == 1
        assert len(b_clusters[1].ions) == 1
        assert b_clusters[1].inclusive_start == 2
        assert b_clusters[1].exclusive_end == 4


class Test_get_y_clusters:
    @staticmethod
    def test_smoke():
        # Arrange
        ions = [
            # y-cluster 1
            DbKmer(
                protein_id=0,
                inclusive_start=1,
                exclusive_end=5,
                charge=2,
                neutral_mass=5.2,
                ion_type=1,
            ),
            DbKmer(
                protein_id=0,
                inclusive_start=3,
                exclusive_end=5,
                charge=2,
                neutral_mass=5.2,
                ion_type=1,
            ),
            # y-cluster 2
            DbKmer(
                protein_id=1,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=1,
            ),
            # b-ions
            DbKmer(
                protein_id=1,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=0,
            ),
            DbKmer(
                protein_id=3,
                inclusive_start=2,
                exclusive_end=4,
                charge=2,
                neutral_mass=5.2,
                ion_type=0,
            ),
        ]

        # Act
        y_clusters = get_y_clusters(ions=ions)

        # Assert
        assert len(y_clusters) == 2

        assert y_clusters[0].protein_id == 0
        assert len(y_clusters[0].ions) == 2
        assert y_clusters[0].inclusive_start == 1
        assert y_clusters[0].exclusive_end == 5

        assert y_clusters[1].protein_id == 1
        assert len(y_clusters[1].ions) == 1
        assert y_clusters[1].inclusive_start == 2
        assert y_clusters[1].exclusive_end == 4
