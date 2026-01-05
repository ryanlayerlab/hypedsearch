import subprocess
import sys
from pathlib import Path

from click.testing import CliRunner

from src.comet_utils import CometPSM
from src.constants import HUMAN_PROTEOME, MOUSE_PROTEOME, RUN_HYPEDSEARCH_SMK
from src.hypedsearch import (
    HybridFormer,
    HybridPSMScorer,
    HypedsearchOutputs,
    HypedsearchRunConfig,
)
from src.hypedsearch_run_analysis import (
    Experiment,
    SpectrumPSMs,
    accept_hybrid_psms_via_neo_fusion,
    cli_create_spectrum_psms,
    find_possible_hybrids_for_seq,
    get_protein_abundance_from_spectra_results,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.peptide_spectrum_comparison import PSM
from src.utils import load_json
from tests.conftest import default_hs_run, default_test_config

FXN4_DIR = Path("tests/data/BMEM_AspN_Fxn4")


class Test_SpectrumCometResults:
    def test_smoke():
        config_path = "results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/3_6_Meoh_RAT_Islet_B35spike/hs.config.json"
        hs_config = HypedsearchRunConfig.from_json(config_path)
        outs = HypedsearchOutputs(hs_config=hs_config)
