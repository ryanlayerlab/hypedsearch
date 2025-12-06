import subprocess
import sys
from pathlib import Path

from click.testing import CliRunner

from src.comet_utils import CometPSM
from src.constants import MOUSE_PROTEOME, RUN_HYPEDSEARCH_SMK
from src.hypedsearch import HybridFormer, HybridPSMScorer, HypedsearchRunConfig
from src.hypedsearch_run_analysis import (
    ExperimentPSMs,
    SpectrumCometResults,
    SpectrumPSMs,
    cli_create_spectrum_psms,
    find_possible_hybrids,
    get_protein_abundance_from_spectra_results,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.peptide_spectrum_comparison import PSM
from src.utils import load_json
from tests.conftest import default_config, default_hs_run

FXN4_DIR = Path("tests/data/BMEM_AspN_Fxn4")


def BMEM_AspN_Fxn4_config(test_data_dir, tmp_path):
    return HypedsearchRunConfig(
        mzml_to_scans={"Fxn": "all"},
        parent_out_dir=tmp_path,
        hybrid_former=HybridFormer(
            kmer_db="n/a",
            fasta=MOUSE_PROTEOME,
        ),
        psm_scorer=HybridPSMScorer(
            comet_params=test_data_dir / "comet.params",
        ),
    )


class Test_get_protein_abundance_from_spectra_results:
    @staticmethod
    def test_smoke(test_data_dir):
        results = SpectrumCometResults.from_dir(
            path=test_data_dir / "BMEM_AspN_Fxn4", by_spectrum=False
        )
        prot_ab = get_protein_abundance_from_spectra_results(
            spectra_results=results,
            q_val_thresh=0.01,
        )


class Test_SpectrumPSMs:
    class Test_from_native_and_hybrid_dirs:
        @staticmethod
        def test_native_and_hybrid_dir_the_same(test_data_dir, tmp_path, snapshot):
            hs_config = default_hs_run(test_data_dir=test_data_dir, out_dir=tmp_path)
            spectra = Mzml(mzml=list(hs_config.mzml_to_scans.keys())[0]).id_to_spectrum
            spectra_psms = SpectrumPSMs.from_native_and_hybrid_dirs(
                native_dir=hs_config.native_run_dir,
                hybrid_dir=hs_config.hybrid_run_dir,
                spectra=spectra,
            )
            for spectrum_psms in spectra_psms:
                assert spectrum_psms.native_target is not None


class Test_cli_create_spectrum_psms:
    @staticmethod
    def test_single_config(test_data_dir):
        # Arrange
        hs_config = test_data_dir / "hs_results/hs_config.json"
        # Act
        runner = CliRunner()
        result = runner.invoke(
            cli_create_spectrum_psms,
            [
                "--config",
                f"{hs_config}",
                "--min_side_len",
                "5",
            ],
        )
        # Assert
        assert result.exit_code == 0
        psms = SpectrumPSMs.load(
            path=SpectrumPSMs.default_save_path(
                out_dir=test_data_dir / "hs_results", min_side_len=5
            )
        )
        assert len(psms) > 0


class Test_find_possible_hybrids:
    @staticmethod
    def test_hybrids_have_right_side_len(tmp_path):
        """ """
        # Arrange
        mzml = "data/251028_RP_Islet_Spikes_Crashout/1_1_Acet_Aspn_Islet_B35spike.mzML"
        db_path = "results/251028_RP_Islet_Spikes_Crashout/inputs/kmers.db"
        scan_num = 1316004
        seq = "GITLNHLKATPIESHQV"
        spectrum = Mzml(mzml=mzml).get_spectrum(scan=scan_num)
        kmer_to_proteins_map = kmer_to_proteins_map = KmerDatabase(
            db_path=db_path
        ).kmer_to_proteins_map.kmer_to_protein_map
        # Act
        hybrids_2 = find_possible_hybrids(
            seq=seq, kmer_to_proteins_map=kmer_to_proteins_map, min_side_len=2
        )
        hybrids_5 = find_possible_hybrids(
            seq=seq, kmer_to_proteins_map=kmer_to_proteins_map, min_side_len=5
        )
        # Assert
        assert len(hybrids_2) > 0
        assert len(hybrids_5) == 0
        assert len(hybrids_5) == 0
        assert len(hybrids_5) == 0
        assert len(hybrids_5) == 0


class Test_ExperimentPSMs:
    class Test_NeoFusion:
        @staticmethod
        def test_smoke():
            # Arrange
            spectra_psms_path = (
                "results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/"
                + "1_1_Acet_Aspn_Islet_B35spike/spectra_psms_minSideLen5.pklz"
            )
            psms = SpectrumPSMs.load(path=spectra_psms_path)
            # Act
            exp = ExperimentPSMs(psms=psms)
            best_iteration, accepted_psms = exp.accept_psm_via_neo_fusion()
