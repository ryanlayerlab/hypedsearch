import json
from dataclasses import asdict

from src.kmer_database import KmerDatabase, KmerToProteinMap, create_kmer_database
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Peptide, UnpositionedProductIon
from src.utils import flatten_list_of_lists, read_new_line_separated_file


class Test_KmerToProteinMap:
    class Test_get_uniq_kmer_to_protein_map:
        @staticmethod
        def test_default():
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
            actual = KmerToProteinMap.get_uniq_kmer_to_protein_map(
                min_k=min_k, max_k=max_k, proteins=proteins
            )

            assert actual == expected

        @staticmethod
        def test_use_name_as_protein_attr():
            proteins = [
                Peptide(seq="ACD", id=0, name="prot 1"),
                Peptide(seq="CDE", id=1, name="prot 2"),
            ]
            min_k, max_k = 1, 3
            expected = {
                "A": ["prot 1"],
                "AC": ["prot 1"],
                "ACD": ["prot 1"],
                "C": ["prot 1", "prot 2"],
                "CD": ["prot 1", "prot 2"],
                "CDE": ["prot 2"],
                "D": ["prot 1", "prot 2"],
                "DE": ["prot 2"],
                "E": ["prot 2"],
            }
            actual = KmerToProteinMap.get_uniq_kmer_to_protein_map(
                min_k=min_k, max_k=max_k, proteins=proteins, protein_attr="name"
            )

            assert actual == expected

        @staticmethod
        def test_kmer_appears_multiple_times_in_one_protein():
            proteins = [
                Peptide(seq="AAA", id=0, name="prot 1"),
            ]
            min_k, max_k = 1, 1
            expected = {"A": ["prot 1"]}
            actual = KmerToProteinMap.get_uniq_kmer_to_protein_map(
                min_k=min_k, max_k=max_k, proteins=proteins, protein_attr="name"
            )
            assert actual == expected

    class Test_create:
        @staticmethod
        def test_from_fasta(test_data_dir):
            fasta_path = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta_path, min_k=1, max_k=3
            )
            assert len(kmer_to_prot_map.kmer_to_protein_map) > 0

        @staticmethod
        def test_from_peptides(test_data_dir):
            peptides = Peptide.from_fasta(
                fasta_path=test_data_dir / "three_proteins.fasta"
            )
            kmer_to_prot_map = KmerToProteinMap.create(
                proteins=peptides, min_k=1, max_k=3
            )
            assert len(kmer_to_prot_map.kmer_to_protein_map) > 0

        @staticmethod
        def test_include_only_proteins_from_file(test_data_dir):
            fasta = (
                test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
            )
            protein_names = test_data_dir / "mouse_data_top_10_proteins.txt"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta,
                min_k=1,
                max_k=3,
                protein_names=protein_names,
                protein_attr="name",
            )
            proteins = set(read_new_line_separated_file(path=protein_names))
            assert proteins == set(
                flatten_list_of_lists(kmer_to_prot_map.kmer_to_protein_map.values())
            )
            # assert len(kmer_to_prot_map.kmer_to_protein_map) > 0
            # # Check that the proteins are from the file
            # for protein in kmer_to_prot_map.kmer_to_protein_map.values():
            #     assert protein[0].startswith("test_")

    class Test_save:
        @staticmethod
        def test_pklz(tmp_path, test_data_dir):
            fasta_path = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta_path, min_k=1, max_k=3
            )
            out_path = tmp_path / "out.pklz"
            kmer_to_prot_map.save(out_path=out_path)
            assert out_path.exists()
            assert len(KmerToProteinMap.load(out_path).kmer_to_protein_map) > 0

        @staticmethod
        def test_json(tmp_path, test_data_dir):
            fasta_path = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta_path, min_k=1, max_k=3
            )
            out_path = tmp_path / "out.json"
            kmer_to_prot_map.save(out_path=out_path)
            assert out_path.exists()
            assert len(KmerToProteinMap.load(out_path).kmer_to_protein_map) > 0


class Test_KmerDatabase:
    class Test_create_db:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path, snapshot, snapshot_dir):
            # Arrange
            fasta_path = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta_path, min_k=1, max_k=25, protein_attr="name"
            )

            # Act
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "three_proteins.db",
                kmer_to_protein_map=kmer_to_prot_map.kmer_to_protein_map,
            )

            # Arrange
            assert kmer_db.db.indices() == [kmer_db.index_name]
            db_rows = kmer_db.db.all_table_rows(kmer_db.table_name)
            assert len(db_rows) == len(kmer_to_prot_map.kmer_to_protein_map)
            # Since we grab the unique kmers in the FASTA, sorting by sequence should
            # produce a unique, reproducible order
            snapshot.snapshot_dir = snapshot_dir
            snapshot_file = "three_proteins.fasta.db.rows"
            db_rows = sorted(db_rows, key=lambda row: row["seq"])
            snapshot.assert_match(
                json.dumps(db_rows, indent=2, sort_keys=True), snapshot_file
            )

    class Test_get_matching_product_ions:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path):
            # Arrange
            peak_mz = 720.3775024414062
            # Create DB
            fasta = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta, min_k=1, max_k=25, protein_attr="name"
            )
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "three_proteins.db",
                kmer_to_protein_map=kmer_to_prot_map.kmer_to_protein_map,
            )

            # Act
            results = kmer_db.get_matching_product_ions(
                query_mz=peak_mz,
                charge=2,
                ppm_tolerance=20,
                ion_type="y",
            )

            # Arrange
            # Check that the results contain at least one expected ion
            truthy_thing = ("AAPAAGSAPAAAEEKK", 2) in {
                (ion.seq, ion.charge) for ion in results
            }
            assert truthy_thing
            # Check that results are product ions
            for ion in results:
                assert isinstance(ion, UnpositionedProductIon)

    class Test_get_peak_ion_matches_for_spectrum:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path, snapshot, snapshot_dir):
            # Arrange
            mzml, scan = test_data_dir / "BMEM_AspN_Fxn4.mzML", 7
            spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
            fasta = (
                test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
            )
            kmer_to_prot_map = KmerToProteinMap.create(
                fasta=fasta,
                min_k=1,
                max_k=25,
                protein_attr="name",
                protein_names=[
                    "sp|P99027|RLA2_MOUSE"
                ],  # using this protein because it's the top Comet PSM for this scan comes from this protein
            )
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "test.db",
                kmer_to_protein_map=kmer_to_prot_map.kmer_to_protein_map,
            )

            # Act
            peak_ion_matches = kmer_db.get_peak_ion_matches_for_spectrum(
                spectrum=spectrum, ppm_tolerance=20
            )

            # Assert
            peak_ion_matches = [asdict(p) for p in peak_ion_matches]
            peak_ion_matches_sorted = sorted(
                peak_ion_matches,
                key=lambda x: (
                    x["peak"]["id"],
                    x["ion"]["seq"],
                    x["ion"]["charge"],
                    x["ion"]["ion_type"],
                ),
            )
            snapshot.snapshot_dir = snapshot_dir
            snapshot_file = "BMEM_AspN_Fxn4_scan7_peak_ion_matches.json"
            snapshot.assert_match(
                json.dumps(peak_ion_matches_sorted, indent=2, sort_keys=True),
                snapshot_file,
            )


class Test_create_db:
    @staticmethod
    def test_fasta_restricted_proteins(test_data_dir, tmp_path, snapshot, snapshot_dir):
        # Arrange
        fasta = test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
        proteins = test_data_dir / "mouse_data_top_10_proteins.txt"
        # Act
        kmer_db = create_kmer_database(
            fasta=fasta,
            kmer_to_proteins_path=tmp_path / "test.pklz",
            db_path=tmp_path / "test.db",
            proteins=proteins,
        )

        # Arrange
        assert kmer_db.db.indices() == [kmer_db.index_name]
        snapshot.snapshot_dir = snapshot_dir
        snapshot_file = f"{fasta.stem}_{proteins.stem}.json"
        # Since we grab the unique kmers in the FASTA, sorting by sequence should
        # produce a unique, reproducible order
        db_rows = kmer_db.db.all_table_rows(kmer_db.table_name)
        db_rows = sorted(db_rows, key=lambda row: row["seq"])
        snapshot.assert_match(
            json.dumps(db_rows, indent=2, sort_keys=True), snapshot_file
        )
