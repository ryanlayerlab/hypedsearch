import subprocess
import sys
from pathlib import Path

from click.testing import CliRunner

from src.comet_utils import CometPSM
from src.constants import HUMAN_PROTEOME, MOUSE_PROTEOME, RUN_HYPEDSEARCH_SMK
from src.hypedsearch import HybridFormer, HybridPSMScorer, HypedsearchRunConfig
from src.hypedsearch_run_analysis import HypedsearchOutputs
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.peptide_spectrum_comparison import PSM
from src.utils import flatten_list_of_lists, load_json
from tests.conftest import default_hs_run, default_test_config

FXN4_DIR = Path("tests/data/BMEM_AspN_Fxn4")


class Test_SpectrumCometResults:
    def test_smoke():
        config_path = "results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/3_6_Meoh_RAT_Islet_B35spike/hs.config.json"
        hs_config = HypedsearchRunConfig.from_json(config_path)
        outs = HypedsearchOutputs(hs_config=hs_config)


class Test_HypedsearchOutputs:
    @staticmethod
    def test_remove_methylation():
        config_path = "results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/3_6_Meoh_RAT_Islet_B35spike/hs.config.json"
        hs_config = HypedsearchRunConfig.from_json(config_path)
        outs = HypedsearchOutputs(hs_config=hs_config, remove_carbamidomethylation=True)
        for hybrid in flatten_list_of_lists(outs.seq_to_hybrids_map.values()):
            assert not hybrid.evidence_of_carbamidomethylation

    @staticmethod
    def test_get_native_and_hybrid_results():
        config_path = "results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/3_6_Meoh_RAT_Islet_B35spike/hs.config.json"
        hs_config = HypedsearchRunConfig.from_json(config_path)
        outs = HypedsearchOutputs(hs_config=hs_config, remove_carbamidomethylation=True)
        spectrum_uid = list(outs._get_hybrid_targets.keys())[0]
        # Assert
        outs.get_native_results_for_spectrum(spectrum_uid=spectrum_uid)
        outs.get_hybrid_results_for_spectrum(spectrum_uid=spectrum_uid)
