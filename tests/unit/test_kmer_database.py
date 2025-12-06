import json
from dataclasses import asdict

from src.kmer_database import DbKmer, KmerDatabase, KmerToProteinsMap
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Fasta, Peptide, UnpositionedProductIon
from src.utils import flatten_list_of_lists, load_json, read_new_line_separated_file


class Test_DbKmer:
    @staticmethod
    def test_check_proteins_are_sorted():
        db_kmer = DbKmer.from_seq_and_proteins(
            seq="ACD", proteins=["protB", "protA", "protB"]
        )
        assert db_kmer.proteins_as_set == {"protA", "protB"}


class Test_KmerToProteinMap:
    class Test_get_uniq_kmer_to_protein_map:
        @staticmethod
        def test_default():
            proteins = [Peptide(seq="ACD", id=0), Peptide(seq="CDE", id=1)]
            min_k, max_k = 1, 3
            expected = {
                "A": {0},
                "AC": {0},
                "ACD": {0},
                "C": {0, 1},
                "CD": {0, 1},
                "CDE": {1},
                "D": {0, 1},
                "DE": {1},
                "E": {1},
            }
            actual = KmerToProteinsMap.get_uniq_kmer_to_protein_map(
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
                "A": {"prot 1"},
                "AC": {"prot 1"},
                "ACD": {"prot 1"},
                "C": {"prot 1", "prot 2"},
                "CD": {"prot 1", "prot 2"},
                "CDE": {"prot 2"},
                "D": {"prot 1", "prot 2"},
                "DE": {"prot 2"},
                "E": {"prot 2"},
            }
            actual = KmerToProteinsMap.get_uniq_kmer_to_protein_map(
                min_k=min_k, max_k=max_k, proteins=proteins, protein_attr="name"
            )

            assert actual == expected

        @staticmethod
        def test_kmer_appears_multiple_times_in_one_protein():
            proteins = [
                Peptide(seq="AAA", id=0, name="prot 1"),
            ]
            min_k, max_k = 1, 1
            expected = {"A": {"prot 1"}}
            actual = KmerToProteinsMap.get_uniq_kmer_to_protein_map(
                min_k=min_k, max_k=max_k, proteins=proteins, protein_attr="name"
            )
            assert actual == expected

    class Test_create:
        @staticmethod
        def test_from_fasta(test_data_dir):
            fasta_path = test_data_dir / "three_proteins.fasta"
            kmer_to_prot_map = KmerToProteinsMap.create(
                fasta=fasta_path, min_k=1, max_k=3
            )
            assert len(kmer_to_prot_map.kmer_to_protein_map) > 0

        @staticmethod
        def test_from_peptides(test_data_dir):
            peptides = Peptide.from_fasta(
                fasta_path=test_data_dir / "three_proteins.fasta"
            )
            kmer_to_prot_map = KmerToProteinsMap.create(
                proteins=peptides, min_k=1, max_k=3
            )
            assert len(kmer_to_prot_map.kmer_to_protein_map) > 0

        @staticmethod
        def test_include_only_proteins_from_file(test_data_dir):
            fasta = (
                test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
            )
            protein_names = test_data_dir / "mouse_data_top_10_proteins.txt"
            kmer_to_prot_map = KmerToProteinsMap.create(
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

    @staticmethod
    def test_save_and_load_pklz(tmp_path, test_data_dir):
        fasta_path = test_data_dir / "three_proteins.fasta"
        kmer_to_prot_map = KmerToProteinsMap.create(fasta=fasta_path, min_k=1, max_k=3)
        out_path = tmp_path / "out.pklz"
        kmer_to_prot_map.save(out_path=out_path)
        loaded_map = KmerToProteinsMap.load(path=out_path)
        assert kmer_to_prot_map == loaded_map

    @staticmethod
    def test_save_and_load_json(tmp_path, test_data_dir):
        fasta_path = test_data_dir / "three_proteins.fasta"
        kmer_to_prot_map = KmerToProteinsMap.create(
            fasta=fasta_path, min_k=1, max_k=3, protein_attr="name"
        )
        out_path = tmp_path / "out.json"
        kmer_to_prot_map.save(out_path=out_path)
        loaded_map = KmerToProteinsMap.load(path=out_path)
        assert kmer_to_prot_map == loaded_map


class Test_KmerDatabase:
    @staticmethod
    def test_get_kmer_to_proteins_map_from_db(test_data_dir):
        # Arrange
        db_path = (
            test_data_dir / "sp-P99027-RLA2_MOUSE_mzml=BMEM_AspN_Fxn4;scan=7_kmer_db.db"
        )
        kmer_db = KmerDatabase(db_path=db_path)
        true_kmer_to_proteins = KmerToProteinsMap.load(
            path=test_data_dir
            / "sp-P99027-RLA2_MOUSE_mzml=BMEM_AspN_Fxn4;scan=7_kmer_to_proteins.json"
        )
        # Act
        computed_kmer_to_proteins = kmer_db.kmer_to_proteins_map
        assert computed_kmer_to_proteins == true_kmer_to_proteins

    class Test_create_db:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path, snapshot, snapshot_dir):
            # Arrange
            fasta_path = test_data_dir / "three_proteins.fasta"
            min_k, max_k = 1, 3
            # Act
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "three_proteins.db",
                proteins=Fasta(path=fasta_path).proteins,
                min_k=min_k,
                max_k=max_k,
            )

            # Assert
            kmer_to_prot_map = KmerToProteinsMap.create(
                fasta=fasta_path, min_k=min_k, max_k=max_k
            )
            assert kmer_db.db.indices() == [kmer_db.index_name]
            db_rows = kmer_db.get_all_rows(as_dicts=True)
            assert len(db_rows) == len(kmer_to_prot_map.kmer_to_protein_map)
            # Since we grab the unique kmers in the FASTA, sorting by sequence should
            # produce a unique, reproducible order

            db_rows = sorted(db_rows, key=lambda row: row["seq"])
            snapshot.snapshot_dir = snapshot_dir
            snapshot_file = "three_proteins.fasta.db.rows"
            snapshot_data = load_json(path=snapshot_dir / snapshot_file)
            for row_idx, row in enumerate(db_rows):
                snapshot_row = snapshot_data[row_idx]
                assert row["seq"] == snapshot_row["seq"]
                assert row["aa_mass"] == snapshot_row["aa_mass"]
                assert (
                    DbKmer(**row).proteins_as_set
                    == DbKmer(**snapshot_row).proteins_as_set
                )

    class Test_get_matching_product_ions:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path):
            # Arrange
            peak_mz = 720.3775024414062
            # Create DB
            fasta = test_data_dir / "three_proteins.fasta"
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "three_proteins.db",
                proteins=Fasta(path=fasta).proteins,
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
            mzml, scan = test_data_dir / "BMEM_AspN_Fxn4/BMEM_AspN_Fxn4.mzML", 7
            spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
            fasta = (
                test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
            )
            protein_names = ["sp|P99027|RLA2_MOUSE"]
            kmer_db = KmerDatabase.create_db(
                db_path=tmp_path / "test.db",
                proteins=Fasta(path=fasta).get_proteins_by_name(names=protein_names),
                min_k=1,
                max_k=25,
            )

            # Act
            peak_ion_matches = kmer_db.get_peak_ion_matches_for_spectrum(
                spectrum=spectrum, ppm_tolerance=20, max_allowed_ion_charge=4
            )

            # Assert
            # Doing a snapshot test because there are 142 peak-ion matches which is too
            # many to individually test
            peak_ion_matches = [p.model_dump() for p in peak_ion_matches]
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
