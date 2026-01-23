from pathlib import Path

from src.constants import HUMAN_PROTEOME
from src.hybrid_simulation import (
    create_new_kmer_db_and_fasta_for_hybrid_simulation,
    make_native_seq_a_hybrid_seq_in_prot,
    prepare_psm_for_hybrid_finding_simulation_study,
    run_hybrid_simulation_on_psm,
)
from src.hypedsearch import (
    HybridFormer,
    HypedsearchRunConfig,
    SpectrumPreprocessor,
    hybrid_run_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta
from src.psm import CometPSM
from tests.conftest import default_test_config


class Test_make_native_seq_a_hybrid_seq_in_prot:
    def test_smoke(mouse_fasta):
        # Arrange
        seq, prot_seq = "AB", "WXABYZ"
        exp_left, exp_right, exp_prot = "A", "B", "WAXYBZ"
        # Act
        left_hy_seq, right_hy_seq, new_prot_seq = make_native_seq_a_hybrid_seq_in_prot(
            native_seq=seq, prot_seq=prot_seq
        )
        # Assert
        assert left_hy_seq == exp_left
        assert right_hy_seq == exp_right
        assert new_prot_seq == exp_prot


def test_create_comet_psm(test_data_dir):
    psms = CometPSM.from_txt(
        txt="results/hs_mouse_samples/native_run/assign-confidence.txt"
    )
    best_psm = max(psms, key=lambda obj: obj.xcorr)
    best_psm.save(path=test_data_dir / "example_comet_psm.json")
    assert 0 == 0


class Test_prepare_psm_for_hybrid_finding_simulation_study:
    @staticmethod
    def test_psm_that_works(mouse_fasta, comet_psm):
        # Arrange
        exp_prot = comet_psm.proteins[0]
        original_prot_seq = Fasta(path=mouse_fasta).protein_name_to_seq_map[exp_prot]
        result = prepare_psm_for_hybrid_finding_simulation_study(
            psm=comet_psm, fasta=mouse_fasta
        )
        assert result.seq in original_prot_seq
        assert result.seq not in result.new_prot_seq
        assert result.left_hy_seq in result.new_prot_seq
        assert result.right_hy_seq in result.new_prot_seq


class Test_create_new_kmer_db_and_fasta_for_hybrid_simulation:
    @staticmethod
    def test_smoke(mouse_fasta, comet_psm, tmp_path):
        # Arrange
        native_to_hybrid_change = prepare_psm_for_hybrid_finding_simulation_study(
            psm=comet_psm, fasta=mouse_fasta
        )
        existing_kmer_db = Path("results/hs_mouse_samples/kmers.db")
        new_kmer_db = tmp_path / "kmers.db"
        new_fasta = tmp_path / "proteins.fasta"
        # Act
        kmer_db = create_new_kmer_db_and_fasta_for_hybrid_simulation(
            native_to_hybrid_change=native_to_hybrid_change,
            fasta=mouse_fasta,
            existing_kmer_db=existing_kmer_db,
            new_kmer_db=new_kmer_db,
            new_fasta=new_fasta,
        )
        kmers_in_db = set(kmer_db.kmer_to_proteins_map.kmer_to_protein_map.keys())
        # Assert
        assert native_to_hybrid_change.prot_name in kmer_db.proteins
        assert native_to_hybrid_change.seq not in kmers_in_db
        assert native_to_hybrid_change.left_hy_seq in kmers_in_db
        assert native_to_hybrid_change.right_hy_seq in kmers_in_db


class Test_run_hybrid_simulation_on_psm:
    @staticmethod
    def test_mouse_psm(mouse_fasta, comet_psm, tmp_path, test_data_dir):
        # Arrange
        spectrum = Spectrum.get_spectrum(
            scan=2376, mzml=test_data_dir / "BMEM_AspN_Fxn4/BMEM_AspN_Fxn4.mzML"
        )
        # Act
        comet_outputs = run_hybrid_simulation_on_psm(
            psm=comet_psm,
            spectrum=spectrum,
            fasta=mouse_fasta,
            out_dir=tmp_path,
            kmer_db=Path("results/hs_mouse_samples/kmers.db"),
            comet_params=Path("results/hybrid_simulation/inputs/crux.comet.params"),
        )
        # Assert
        hpsm = [
            psm for psm in CometPSM.from_txt(txt=comet_outputs.target) if psm.num == 1
        ][0]
        assert hpsm.seq == comet_psm.seq

    @staticmethod
    def test_human_psm(tmp_path):
        # Arrange
        q_thresh = 0.01
        psms = CometPSM.from_txt(
            txt="results/1_1_Acet_Aspn_Islet_B35spike/native_run/assign-confidence.target.txt"
        )
        psms = [psm for psm in psms if psm.q_value <= q_thresh]
        psm = psms[0]
        mzml = Mzml(
            path="data/251028_RP_Islet_Spikes_Crashout/1_1_Acet_Aspn_Islet_B35spike.mzML"
        )
        kmer_db_path = Path("results/hybrid_simulation/inputs/kmers.db")
        # Act
        comet_outputs = run_hybrid_simulation_on_psm(
            psm=psm,
            spectrum=mzml.get_spectrum(scan=psm.scan),
            fasta=HUMAN_PROTEOME,
            out_dir=tmp_path,
            kmer_db=kmer_db_path,
            comet_params=Path("results/hybrid_simulation/inputs/crux.comet.params"),
        )
        # Assert
        hpsm = [
            psm for psm in CometPSM.from_txt(txt=comet_outputs.target) if psm.num == 1
        ][0]
        assert hpsm.seq == psm.seq
