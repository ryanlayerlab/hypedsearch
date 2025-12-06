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
from src.crux import CometConfig
from src.constants import TARGET

smk_config = SimpleNamespace(**config)
comet_config = CometConfig.from_json(path=smk_config.comet_config)
expected_outputs = comet_config.expected_outputs(psm_type=TARGET)
missing_outputs = comet_config.missing_outputs(psm_type=TARGET)

print(f"Given the CometConfig, there are {len(expected_outputs)} of which {len(missing_outputs)} are missing.")

rule all:
    input:
        missing_outputs

rule run_comet:
    input: 
        mzml = lambda wildcards: next(
            Path(mzml) for mzml in comet_config.mzml_to_scans.keys() if Mzml.get_mzml_name(mzml) == wildcards.sample
        ),
    output:
        target = "{out_dir}/{sample}.comet.{scan}-{scan}.target.txt"
    singularity: 
        "docker://airikjohnson/hypedsearch:latest"
    shell:
        """ 
        python -m src.crux run-comet \
            --mzml {input.mzml} \
            --config {smk_config.comet_config} \
            --scan_min {wildcards.scan} \
            --scan_max {wildcards.scan}
        """