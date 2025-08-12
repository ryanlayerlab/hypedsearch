import json
from dataclasses import asdict

import pytest

from src.constants import B_ION_TYPE, Y_ION_TYPE
from src.hybrids_via_clusters import (
    Cluster,
    PositionedProductIon,
    SpectrumClusters,
    form_extended_clusters_for_spectrum,
    form_hybrids_from_clusters,
    form_spectrum_hybrids_via_clustering,
)
from src.kmer_database import create_db
from src.mass_spectra import Spectrum
from src.peptides_and_ions import UnpositionedProductIon, compute_peptide_precursor_mz
from src.utils import mass_difference_in_ppm

# def test_testing_shit(test_data_dir, tmp_path):
# Arrange
# mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
# spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
# db_path = tmp_path / "test.db"
# fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
# kmer_db = create_db(
#     kmer_to_protein_path=tmp_path / "test.pklz",
#     fasta=test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta",
#     proteins=[
#         "sp|P99027|RLA2_MOUSE"
#     ],  # using this protein because it's the top Comet PSM for this scan comes from this protein,
#     db_path=db_path,
# )


class Test_PositionedProductIon:
    @staticmethod
    def test_from_unpositioned_product_ion():
        # Arrange
        ion = UnpositionedProductIon(
            seq="ABC", charge=2, ion_type=B_ION_TYPE, proteins=["prot1", "prot2"]
        )
        prot_name_to_seq = {"prot1": "XXABCXXABCXX", "prot2": "ABCXX"}
        # Act
        result = PositionedProductIon.from_unpositioned_product_ion(
            unpositioned_product_ion=ion, protein_name_to_seq_map=prot_name_to_seq
        )
        # Assert
        expected = [
            PositionedProductIon(
                seq="ABC",
                charge=2,
                ion_type=B_ION_TYPE,
                protein="prot1",
                inclusive_start=2,
                exclusive_end=5,
            ),
            PositionedProductIon(
                seq="ABC",
                charge=2,
                ion_type=B_ION_TYPE,
                protein="prot1",
                inclusive_start=7,
                exclusive_end=10,
            ),
            PositionedProductIon(
                seq="ABC",
                charge=2,
                ion_type=B_ION_TYPE,
                protein="prot2",
                inclusive_start=0,
                exclusive_end=3,
            ),
        ]
        assert result == expected


class Test_SpectrumClusters:
    @staticmethod
    def test_get_clusters():
        # Arrange
        ions = [
            PositionedProductIon(
                seq="ABC",
                charge=2,
                ion_type=B_ION_TYPE,
                protein="prot1",
                inclusive_start=2,
                exclusive_end=5,
            ),
            PositionedProductIon(
                seq="ABCD",
                charge=1,
                ion_type=B_ION_TYPE,
                protein="prot1",
                inclusive_start=2,
                exclusive_end=6,
            ),
            PositionedProductIon(
                seq="AB",
                charge=1,
                ion_type=Y_ION_TYPE,
                protein="prot2",
                inclusive_start=1,
                exclusive_end=3,
            ),
            PositionedProductIon(
                seq="AB",
                charge=1,
                ion_type=B_ION_TYPE,
                protein="prot2",
                inclusive_start=1,
                exclusive_end=3,
            ),
        ]
        # Act
        b_clusters = SpectrumClusters._get_clusters(
            positioned_ions=ions, ion_type=B_ION_TYPE
        )
        y_clusters = SpectrumClusters._get_clusters(
            positioned_ions=ions, ion_type=Y_ION_TYPE
        )
        # Assert
        expected_b_clusters = [
            Cluster(ions=[ions[0], ions[1]]),
            Cluster(ions=[ions[3]]),
        ]
        assert b_clusters == expected_b_clusters
        expected_y_clusters = [
            Cluster(ions=[ions[2]]),
        ]
        assert y_clusters == expected_y_clusters


class Test_Cluster:
    class Test_post_init:
        @staticmethod
        def test_successful_creation():
            # Act
            ions = [
                PositionedProductIon(
                    seq="ABC",
                    charge=2,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
                PositionedProductIon(
                    seq="ABCD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=6,
                ),
            ]
            # Act
            cluster = Cluster(ions=ions)
            # Assert
            assert cluster.ion_type == B_ION_TYPE
            assert cluster.inclusive_start == 2
            assert cluster.exclusive_end == 6

        @staticmethod
        def test_different_ion_types():
            # Act
            ions = [
                PositionedProductIon(
                    seq="ABC",
                    charge=2,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
                PositionedProductIon(
                    seq="ABCD",
                    charge=1,
                    ion_type=Y_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=6,
                ),
            ]
            # Act & Assert
            with pytest.raises(AssertionError, match="same ion type"):
                Cluster(ions=ions)

        @staticmethod
        def test_different_proteins():
            # Act
            ions = [
                PositionedProductIon(
                    seq="ABC",
                    charge=2,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
                PositionedProductIon(
                    seq="ABCD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot2",
                    inclusive_start=2,
                    exclusive_end=6,
                ),
            ]
            # Act & Assert
            with pytest.raises(AssertionError, match="same protein"):
                Cluster(ions=ions)

        @staticmethod
        def test_b_ion_different_starting_points():
            # Act
            ions = [
                PositionedProductIon(
                    seq="ABC",
                    charge=2,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
                PositionedProductIon(
                    seq="ABCD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=3,
                    exclusive_end=6,
                ),
            ]
            # Act & Assert
            with pytest.raises(AssertionError, match="same start position"):
                Cluster(ions=ions)

        @staticmethod
        def test_y_ion_different_end_points():
            # Act
            ions = [
                PositionedProductIon(
                    seq="ABC",
                    charge=2,
                    ion_type=Y_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
                PositionedProductIon(
                    seq="ABCD",
                    charge=1,
                    ion_type=Y_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=6,
                ),
            ]
            # Act & Assert
            with pytest.raises(AssertionError, match="same end position"):
                Cluster(ions=ions)

    class Test_set_extended_seq:
        @staticmethod
        def test_smoke():
            ions = [
                PositionedProductIon(
                    seq="AD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=4,
                ),
                PositionedProductIon(
                    seq="ADN",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
            ]
            cluster = Cluster(ions=ions)
            # Act
            # m/z at charge=2: (ADN, ~160), (ADNE, ~224), (ADNEE, ~289), (ADNEEE, ~353)
            cluster.set_extended_seq(
                protein_seq="EEADNEEE",
                precursor_mz_ppm_tol=20,
                precursor_charge=2,
                precursor_mz=300,
            )
            # Assert
            assert cluster.extended_seq == "ADNEE"

    class Test_get_seqs_for_hybrids:
        @staticmethod
        def test_b_ions():
            ions = [
                PositionedProductIon(
                    seq="AD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=4,
                ),
                PositionedProductIon(
                    seq="NAD",
                    charge=1,
                    ion_type=B_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=5,
                ),
            ]
            cluster = Cluster(ions=ions)
            cluster.extended_seq = "ADNEE"
            # Act
            seqs = cluster.get_seqs_for_hybrids()
            # Assert
            assert seqs == ["AD", "ADN", "ADNE", "ADNEE"]

        @staticmethod
        def test_y_ions():
            ions = [
                PositionedProductIon(
                    seq="AD",
                    charge=1,
                    ion_type=Y_ION_TYPE,
                    protein="prot1",
                    inclusive_start=2,
                    exclusive_end=4,
                ),
                PositionedProductIon(
                    seq="NAD",
                    charge=1,
                    ion_type=Y_ION_TYPE,
                    protein="prot1",
                    inclusive_start=1,
                    exclusive_end=4,
                ),
            ]
            cluster = Cluster(ions=ions)
            cluster.extended_seq = "EENAD"
            # Act
            seqs = cluster.get_seqs_for_hybrids()
            # Assert
            assert seqs == ["EENAD", "ENAD", "NAD", "AD"]


class Test_form_extended_clusters_for_spectrum:
    @staticmethod
    def test_smoke(test_data_dir, tmp_path, snapshot, snapshot_dir):
        # Arrange
        mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
        fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        db_path = tmp_path / "test.db"
        create_db(
            kmer_to_protein_path=tmp_path / "test.pklz",
            fasta=fasta,
            proteins=[
                "sp|P99027|RLA2_MOUSE"
            ],  # using this protein because the top Comet PSM for this scan comes from this protein,
            db_path=db_path,
        )
        spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)

        # Act
        clusters = form_extended_clusters_for_spectrum(
            db_path=db_path, spectrum=spectrum, fasta=fasta
        )
        # Assert
        # Remove some information so comparison is more straightforward
        clusters_for_json = {B_ION_TYPE: [], Y_ION_TYPE: []}
        keys_to_include = (
            "protein",
            "inclusive_start",
            "exclusive_end",
            "extended_seq",
        )
        for cluster in clusters.b_clusters:
            d = asdict(cluster)
            clusters_for_json[B_ION_TYPE].append(
                {key: d[key] for key in keys_to_include}
            )
        for cluster in clusters.y_clusters:
            d = asdict(cluster)
            clusters_for_json[Y_ION_TYPE].append(
                {key: d[key] for key in keys_to_include}
            )
        clusters_for_json[B_ION_TYPE] = sorted(
            clusters_for_json[B_ION_TYPE],
            key=lambda x: (x["protein"], x["inclusive_start"], x["exclusive_end"]),
        )
        clusters_for_json[Y_ION_TYPE] = sorted(
            clusters_for_json[Y_ION_TYPE],
            key=lambda x: (x["protein"], x["inclusive_start"], x["exclusive_end"]),
        )
        snapshot.snapshot_dir = snapshot_dir
        snapshot_file = "BMEM_AspN_Fxn4_scan7_clusters.json"
        snapshot.assert_match(
            json.dumps(clusters_for_json, indent=2, sort_keys=True),
            snapshot_file,
        )


class Test_form_hybrids_from_clusters:
    @staticmethod
    def test_smoke(test_data_dir, tmp_path, snapshot, snapshot_dir):
        # Arrange
        mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
        fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        db_path = tmp_path / "test.db"
        create_db(
            kmer_to_protein_path=tmp_path / "test.pklz",
            fasta=fasta,
            proteins=test_data_dir
            / "mouse_data_top_10_proteins.txt",  # using this protein because the top Comet PSM for this scan comes from this protein,
            db_path=db_path,
        )
        spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
        clusters = form_extended_clusters_for_spectrum(
            db_path=db_path, spectrum=spectrum, fasta=fasta
        )
        # Act
        hybrids = form_hybrids_from_clusters(
            b_clusters=clusters.b_clusters,
            y_clusters=clusters.y_clusters,
            precursor_mz=spectrum.precursor_mz,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz_ppm_tol=20,
            scan=spectrum.scan,
            sample=spectrum.sample,
        )
        # Assert
        # Remove some information so comparison is more straightforward
        hybrids_for_json = []
        for hybrid in hybrids:
            # hybrids_for_json.append(f"{hybrid.left_seq}{hybrid.right_seq}")
            hybrids_for_json.append(
                f"{hybrid.left_seq}{hybrid.right_seq}|{'|'.join(sorted(hybrid.left_proteins))}|{'|'.join(sorted(hybrid.right_proteins))}"
            )
        hybrids_for_json = sorted(hybrids_for_json)
        snapshot.snapshot_dir = snapshot_dir
        snapshot_file = "BMEM_AspN_Fxn4_scan7_hybrids.json"
        snapshot.assert_match(
            json.dumps(hybrids_for_json, indent=2, sort_keys=True),
            snapshot_file,
        )


class Test_form_spectrum_hybrids_via_clustering:
    @staticmethod
    def test_smoke(test_data_dir, tmp_path):
        # Arrange
        mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
        spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
        fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        db_path = tmp_path / "test.db"
        create_db(
            kmer_to_protein_path=tmp_path / "test.json",
            fasta=fasta,
            proteins=test_data_dir
            / "mouse_data_top_10_proteins.txt",  # using this protein because the top Comet PSM for this scan comes from this protein,
            db_path=db_path,
        )
        ppm_tol = 20
        # Act
        seq_to_hybrids = form_spectrum_hybrids_via_clustering(
            database=db_path,
            spectrum=spectrum,
            fasta=fasta,
            precursor_mz_ppm_tol=ppm_tol,
        )
        # Assert
        for seq, hybrids in seq_to_hybrids.items():
            assert (
                mass_difference_in_ppm(
                    mass1=spectrum.precursor_mz,
                    mass2=compute_peptide_precursor_mz(
                        seq=seq, charge=spectrum.precursor_charge
                    ),
                )
                <= ppm_tol
            )


def test_stuff(test_data_dir, tmp_path):
    mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
    spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
    fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
    db_path = tmp_path / "test.db"
    create_db(
        kmer_to_protein_path=tmp_path / "test.json",
        fasta=fasta,
        proteins=test_data_dir
        / "mouse_data_top_10_proteins.txt",  # using this protein because the top Comet PSM for this scan comes from this protein,
        db_path=db_path,
    )
    ppm_tol = 20
    # Act
    seq_to_hybrids = form_spectrum_hybrids_via_clustering(
        database=db_path,
        spectrum=spectrum,
        fasta=fasta,
        precursor_mz_ppm_tol=ppm_tol,
    )
    pass
