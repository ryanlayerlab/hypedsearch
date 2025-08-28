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


comet_config = CometConfig(**config)
config = SimpleNamespace(**config)
expected_outputs = comet_config.expected_outputs(
    out_dir=config.out_dir, psm_type=TARGET
)
print(f"There are {len(expected_outputs)} expected outputs.")

rule all:
    input:
        expected_outputs

rule run_comet:
    input: 
        mzml = lambda wildcards: next(
            Path(mzml) for mzml in config.mzml_to_scans.keys() if Mzml.get_mzml_name(mzml) == wildcards.sample
        ),
    output:
        target = "{out_dir}/{sample}.comet.{scan}-{scan}.target.txt"
    singularity: 
        "docker://airikjohnson/hypedsearch:latest"
    script:
        "scripts/run_comet.py"
    # run:
    #     from src.crux import Crux
    #     from src.utils import setup_logger
    #     log_file = Path(f"logs/{Mzml.get_mzml_name(input.mzml)}_{wildcards.scan}.log")
    #     logger = setup_logger(str(log_file))
    #     logger.info(f"Logging to {log_file}")
    #     Crux.run_comet(
    #         mzml=input.mzml,
    #         fasta=config.fasta,
    #         crux_comet_params=config.crux_comet_params,
    #         decoy_search=config.decoy_search,
    #         out_dir=config.out_dir,
    #         file_root=wildcards.sample,
    #         scan_min=int(wildcards.scan),
    #         scan_max=int(wildcards.scan)
    #     )
    #     log_file.unlink()

