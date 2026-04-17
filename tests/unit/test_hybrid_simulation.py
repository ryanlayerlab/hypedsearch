import random
from pathlib import Path

from src.constants import MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.hybrid_simulation import (
    RANDOM,
    HybridSimulationExperiment,
    create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots,
    cut_seq_into_hybrid,
    run_hybrid_simulation_on_spectrum,
    validate_hybrid_for_hybrid_simulation,
)
from src.hypedsearch import HybridRunParams
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.peptides_and_ions import Fasta, Peptide
from src.psm import CometPSM
from src.utils import load_json, to_json


class Test_create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots:
    @staticmethod
    def test_smoke(tmp_path, test_data_dir):
        # Arrange
        left_seq, right_seq = "EVEDPQVEQLELGG", "SPGDLQTLALEVARQ"
        prot_containing_seq = "sp|P01325|INS1_MOUSE"
        seq = left_seq + right_seq
        db = KmerDatabase(db_path=test_data_dir / "mouse_top_10_proteins.db")
        fasta = Fasta(path=MOUSE_PROTEOME)
        new_kmer_db_path, new_fasta_path = (
            tmp_path / "new_kmers.db",
            tmp_path / "new_fasta.fasta",
        )
        # Act
        create_new_kmer_db_and_fasta_for_simulation_with_hybrid_as_new_prots(
            existing_kmer_db=db,
            fasta=fasta,
            left_aa_seq=left_seq,
            right_aa_seq=right_seq,
            hybridized_kmer_db=new_kmer_db_path,
            hybridized_fasta=new_fasta_path,
        )
        # Assert
        new_kmer_db = KmerDatabase(db_path=new_kmer_db_path)
        new_fasta = Fasta(path=new_fasta_path)
        # k-mer DB assertions
        assert left_seq in new_kmer_db.kmers
        assert right_seq in new_kmer_db.kmers
        assert seq not in new_kmer_db.kmers
        # FASTA assertions
        assert prot_containing_seq not in new_fasta.protein_name_to_seq_map
        assert f"{prot_containing_seq}|hybridized" in new_fasta.protein_name_to_seq_map
        assert (
            seq
            not in new_fasta.protein_name_to_seq_map[
                f"{prot_containing_seq}|hybridized"
            ]
        )


class Test_cut_seq_into_hybrid:
    @staticmethod
    def test_cut_in_half():
        assert ("AB", "CD") == cut_seq_into_hybrid(seq="ABCD", min_side_len=2)

    @staticmethod
    def test_cut_randomly():
        left_seq, right_seq = cut_seq_into_hybrid(
            seq="ABCDEF", min_side_len=2, method=RANDOM
        )
        assert left_seq + right_seq == "ABCDEF"
        assert len(left_seq) >= 2
        assert len(right_seq) >= 2

    @staticmethod
    def test_cut_at_proportion():
        left_seq, right_seq = cut_seq_into_hybrid(
            seq="ABCDE",
            min_side_len=2,
            method=0.4,
        )
        assert left_seq + right_seq == "ABCDE"

    @staticmethod
    def test_cut_at_idx():
        seq = "ABCDE"
        left_seq, right_seq = cut_seq_into_hybrid(
            seq="ABCDE",
            min_side_len=2,
            method=2,
        )
        assert left_seq == seq[:2]
        assert right_seq == seq[2:]


class Test_validate_hybrid_for_hybrid_simulation:
    @staticmethod
    def test_appears_in_one_protein_once(tmp_path):
        Fasta.write_fasta(
            peptides=[
                Peptide(seq="AXXXY", name="prot1"),
                Peptide(seq="ABC", name="prot2"),
            ],
            path=tmp_path / "proteome.fasta",
        )
        assert validate_hybrid_for_hybrid_simulation(
            left_seq="X",
            right_seq="XX",
            fasta=Fasta(path=tmp_path / "proteome.fasta"),
            min_side_len=1,
        )

    @staticmethod
    def test_appears_in_two_proteins(tmp_path):
        Fasta.write_fasta(
            peptides=[
                Peptide(seq="AXXXY", name="prot1"),
                Peptide(seq="ABXXXC", name="prot2"),
            ],
            path=tmp_path / "proteome.fasta",
        )
        assert not validate_hybrid_for_hybrid_simulation(
            left_seq="X",
            right_seq="XX",
            fasta=Fasta(path=tmp_path / "proteome.fasta"),
            min_side_len=1,
        )

    @staticmethod
    def test_appears_in_one_protein_twice(tmp_path):
        Fasta.write_fasta(
            peptides=[
                Peptide(seq="AXXXYYXXXCD", name="prot1"),
                Peptide(seq="ABC", name="prot2"),
            ],
            path=tmp_path / "proteome.fasta",
        )
        assert not validate_hybrid_for_hybrid_simulation(
            left_seq="X",
            right_seq="XX",
            fasta=Fasta(path=tmp_path / "proteome.fasta"),
            min_side_len=1,
        )


class Test_run_hybrid_simulation_on_spectrum:
    @staticmethod
    def test_smoke(test_data_dir, mouse_mzml_path, tmp_path):
        ndir = tmp_path / "natives"
        ndir.mkdir(parents=True, exist_ok=True)
        hdir = tmp_path / "hybrids"
        mzml = Mzml(path=mouse_mzml_path)
        hdir.mkdir(parents=True, exist_ok=True)
        run_hybrid_simulation_on_spectrum(
            scan=7,
            hybrid_run_params=HybridRunParams.load(
                path=test_data_dir / "hybrid_run_params.json"
            ),
            mzml=mzml,
            left_seq="SAAPAAGS",
            right_seq="APAAAEEKK",
            native_out_dir=ndir,
            hybrid_out_dir=hdir,
            crux_path=MAC_CRUX_EXECUTABLE,
        )
        txt = ndir / f"{mzml.name}.comet.7-7.target.txt"
        assert txt.exists()
        assert len(CometPSM.from_txt(txt=txt)) > 0

        txt = hdir / f"{mzml.name}.comet.7-7.target.txt"
        assert txt.exists()
        h_psms = CometPSM.from_txt(txt=txt)
        assert len(h_psms) > 0
        assert h_psms[0].seq == "SAAPAAGSAPAAAEEKK"


class Test_HybridSimulationExperiment:
    class Test_create_experiment:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir, mouse_mzml_path):
            psms = CometPSM.from_txt(
                txt=test_data_dir / "BMEM_AspN_Fxn4.assign-confidence.txt"
            )
            psms = [psm for psm in psms if psm.q_value < 0.01]
            exp = HybridSimulationExperiment.create_experiment(
                hybrid_run_params=test_data_dir / "hybrid_run_params.json",
                mzml=mouse_mzml_path,
                psms=psms,
                parent_out_dir=tmp_path,
            )
            pass

    class Test_run_experiment_in_parallel:
        @staticmethod
        def test_smoke(tmp_path, mouse_mzml_path, test_data_dir):
            # Arrange
            mzml = Mzml(path=mouse_mzml_path)
            exp = HybridSimulationExperiment(
                hybrid_run_params=HybridRunParams.load(
                    path=test_data_dir / "hybrid_run_params.json"
                ),
                mzml=mzml,
                scan_to_left_right_seq={
                    7: ("SAAPAAGS", "APAAAEEKK"),
                    10: ("SAAP", "AAGSAPAAAEEKK"),
                },
                parent_out_dir=tmp_path,
            )
            # Act
            exp.run_experiment_in_parallel(n_cores=2, crux_path=MAC_CRUX_EXECUTABLE)
            # Assert
            # Native results
            assert len(list(exp.native_out_dir.glob("*.txt"))) == 2
            for txt in exp.native_out_dir.glob("*.txt"):
                assert len(CometPSM.from_txt(txt)) > 0
            # Hybrid results
            assert len(list(exp.hybrid_out_dir.glob("*.txt"))) == 2
            for txt in exp.hybrid_out_dir.glob("*.txt"):
                assert len(CometPSM.from_txt(txt)) > 0


class Test_april_15_2026:
    @staticmethod
    def test_smoke(tmp_path):
        exp = HybridSimulationExperiment.load(
            "results/04-14-26-BMEM_AspN_Fxn4-hybrid-sim-half/hybrid.simulation.config.json"
        )
        # Hybrid sequence is DRVISLSGEHLGRILTGSSEPEAAP
        scan, left_seq, right_seq = 1555, "DRVISLSGEHS", "IIGRTMVVHEKQ"
        native_dir, hybrid_dir = tmp_path / "native", tmp_path / "hybrid"
        native_dir.mkdir(exist_ok=True, parents=True)
        hybrid_dir.mkdir(exist_ok=True, parents=True)
        run_hybrid_simulation_on_spectrum(
            scan=scan,
            left_seq=left_seq,
            right_seq=right_seq,
            mzml=exp.mzml,
            hybrid_run_params=exp.hybrid_run_params,
            native_out_dir=native_dir,
            hybrid_out_dir=hybrid_dir,
            crux_path=MAC_CRUX_EXECUTABLE,
        )
        assert (
            CometPSM.from_txt(list(hybrid_dir.glob("*"))[0])[0].seq
            == left_seq + right_seq
        )
