import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner
from pydantic import ValidationError

from src.constants import (
    HUMAN_PROTEOME,
    MAC_CRUX_EXECUTABLE,
    MOUSE_PROTEOME,
    NATIVE,
    RUN_HYPEDSEARCH_SMK,
    TARGET,
)
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import (
    HybridFormer,
    HybridPSMScorer,
    HypedsearchRunConfig,
    SpectrumPreprocessor,
    SpectrumSelector,
    TrueHybrid,
    cli_run_hypedsearch,
    find_possible_hybrids_for_seq,
    hybrid_run_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta
from src.psm import CometPSM
from src.utils import flatten_list_of_lists, load_json, mass_difference_in_ppm


class Test_HypedsearchRunConfig:
    @staticmethod
    def test_native_run(tmp_path, test_hs_config_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        outputs = config.native_comet_run(
            # dry_run=True
        )
        assert len(outputs) == len(config.mzml_names)
        assert outputs[0].target.exists()
        assert outputs[0].decoy.exists()
        assert len(CometPSM.from_txt(txt=outputs[0].target)) > 0

    @staticmethod
    def test_hybrid_run(tmp_path, test_hs_config_path, mouse_mzml_path, mouse_spectrum):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        seq_to_hybrids, comet_outputs = config.hybrid_run_on_spectrum(
            spectrum=mouse_spectrum
        )
        assert comet_outputs.target.exists()
        assert len(CometPSM.from_txt(txt=comet_outputs.target)) > 0
