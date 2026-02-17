import json
from dataclasses import asdict

import pytest

from src.constants import B_ION_TYPE, MOUSE_PROTEOME, Y_ION_TYPE
from src.hybrids_via_clusters import (
    Cluster,
    HybridPeptide,
    PositionedProductIon,
    SpectrumClusters,
    form_extended_clusters_for_spectrum,
    form_spectrum_hybrids_via_clustering,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Fasta, Fasta2MFMIndex, UnpositionedProductIon


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
            seqs = cluster.get_seqs_for_hybrids(min_side_len=3)
            # Assert
            assert seqs == ["ADN", "ADNE", "ADNEE"]

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
            seqs = cluster.get_seqs_for_hybrids(min_side_len=3)
            # Assert
            assert seqs == ["EENAD", "ENAD", "NAD"]


class Test_form_extended_clusters_for_spectrum:
    @staticmethod
    def test_smoke(snapshot, snapshot_dir, test_data_dir):
        # Arrange
        kmer_db = KmerDatabase(db_path=test_data_dir / "mouse_samples.kmers.db")
        spectrum = Spectrum.get_spectrum(
            mzml=test_data_dir / "mouse_BMEM_AspN_Fxn4.mzML", scan=4538
        )
        min_cluster_len, min_cluster_support = 3, 3
        # Act
        clusters = form_extended_clusters_for_spectrum(
            kmer_db=kmer_db,
            spectrum=spectrum,
            protein_name_to_seq_map={
                prot: Fasta(path=MOUSE_PROTEOME).protein_name_to_seq_map[prot]
                for prot in kmer_db.proteins
            },
            min_cluster_len=min_cluster_len,
            min_cluster_support=min_cluster_support,
            peak_to_ion_ppm_tol=20,
            precursor_mz_ppm_tol=20,
            max_allowed_ion_charge=4,
        )
        # Assert
        for cluster in clusters.b_clusters + clusters.y_clusters:
            assert cluster.support >= min_cluster_support
            assert cluster.length >= min_cluster_len
        snapshot.snapshot_dir = snapshot_dir
        snapshot.assert_match(
            json.dumps(asdict(clusters), indent=2),
            f"{spectrum.uid}.clusters.json",
        )


class Test_form_spectrum_hybrids_via_clustering:
    @staticmethod
    def test_smoke(test_data_dir, snapshot_dir, snapshot):
        # Arrange
        spectrum = Spectrum.get_spectrum(
            mzml=test_data_dir / "mouse_BMEM_AspN_Fxn4.mzML", scan=4538
        )
        seq_to_hybrids = form_spectrum_hybrids_via_clustering(
            spectrum=spectrum,
            kmer_db=KmerDatabase(db_path=test_data_dir / "mouse_samples.kmers.db"),
            fasta=Fasta(path=MOUSE_PROTEOME),
            fasta_fm_index=Fasta2MFMIndex.load(
                path="fastas/SwissProt.TAW_mouse_w_NOD_IAPP.mfmindex"
            ),
            precursor_mz_ppm_tol=20,
            peak_to_ion_ppm_tol=20,
            min_side_len=3,
            min_cluster_support=3,
            max_allowed_ion_charge=4,
            remove_carbamidomethylated_hybrids=True,
        )
        # Assert
        # Normalize + sort before snapshot
        seq_to_hybrids = {
            seq: [
                hy.model_dump(mode="json")
                for hy in sorted(
                    hybrids, key=lambda hy: (len(hy.left_seq), hy.left_seq)
                )
            ]
            for seq, hybrids in seq_to_hybrids.items()
        }
        snapshot.snapshot_dir = snapshot_dir
        snapshot.assert_match(
            json.dumps(seq_to_hybrids, indent=2, sort_keys=True),
            f"{spectrum.uid}.seq2hybrids.json",
        )


class Test_HybridPeptide:
    @staticmethod
    def test_seq_to_hybrids_map_to_peptides_for_fasta():
        # Arrange
        h1 = HybridPeptide(
            left_seq="AB",
            right_seq="C",
            left_proteins=set(["sp|10|INS2_MOUSE", "sp|11|GLUC"]),
            right_proteins=set(["sp|20|INS1_MOUSE"]),
        )
        h2 = HybridPeptide(
            left_seq="A",
            right_seq="BCD",
            left_proteins=set(["sp|30|PROT1"]),
            right_proteins=set(["sp|40|PROT2", "sp|50|PROT3"]),
        )
        h3 = HybridPeptide(
            left_seq="AB",
            right_seq="CD",
            left_proteins=set(["sp|60|PROT4"]),
            right_proteins=set(["sp|70|PROT5"]),
        )
        seq_to_hybrids = {
            "ABC": [h1],
            "ABCD": [h2, h3],
        }
        # Act
        peptides = HybridPeptide.seq_to_hybrids_map_to_peptides(
            seq_to_hybrids=seq_to_hybrids
        )
        assert len(peptides) == 2
        # "ABC" hybrid
        assert peptides[0].seq == "ABC"
        hybrids = HybridPeptide.parse_hybrid_fasta_name(name=peptides[0].name)
        assert len(hybrids) == 1
        assert hybrids[0] == h1
        # "ABCD" hybrids
        assert peptides[1].seq == "ABCD"
        hybrids = HybridPeptide.parse_hybrid_fasta_name(name=peptides[1].name)
        assert len(hybrids) == 2
        assert (hybrids[0] == h2) or (hybrids[1] == h2)
        assert (hybrids[0] == h3) or (hybrids[1] == h3)
