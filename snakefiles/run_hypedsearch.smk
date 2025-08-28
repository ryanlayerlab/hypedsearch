# Expected items in config file:
# 1. mzmls: list of paths to mzML files on which Comet will be run
# 2. out_dir: directory to store the output files
# 3. fasta: path to the native FASTA file
# 4. crux_comet_params: path to the `crux comet` parameter file
# 5. database: path to the k-mer database file from which we'll form hybrids from
# Here's an example config:
# config["mzml_to_scans"] = {
#     "tests/data/mouse_spectra.mzML": [1, 3, 7, 10],
#     "tests/data/mouse_BMEM_AspN_Fxn4.mzML": [2, 4]
# }  
# config["out_dir"] = "tmp/test_hs_run"
# config["fasta"] = "tests/data/mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"
# config["crux_comet_params"] = "tests/data/crux.comet.params"
# config["database"] = "tests/data/mouse_top_10_proteins.db"

import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from typing import List
from types import SimpleNamespace
from src.mass_spectra import Mzml
from src.constants import TARGET
from src.hypedsearch import HybridRunConfig
from src.crux import get_expected_comet_outputs

hs_config = HybridRunConfig(**config)
config = SimpleNamespace(**config)
expected_outputs = hs_config.expected_outputs
print(f"There are {len(expected_outputs)} expected output files")

rule all:
    input:
        expected_outputs

rule run_hypedsearch:
    input: 
        mzml = lambda wildcards: next(
            Path(mzml) for mzml in config.mzml_to_scans.keys() if Mzml.get_mzml_name(mzml) == wildcards.sample
        ),
    output:
        target = f"{config.out_dir}/{{sample}}.comet.{{scan}}-{{scan}}.target.txt"
    benchmark:
        "logs/hypedsearch/{sample}.{scan}.log"
    singularity: 
        "docker://airikjohnson/hypedsearch:latest"
    script:
        "scripts/run_hypedsearch.py"