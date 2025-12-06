import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from typing import List
from types import SimpleNamespace
from src.mass_spectra import Mzml
from src.hypedsearch import HypedsearchRunConfig
from src.utils import setup_logger

# logger = setup_logger()
config = SimpleNamespace(**config)
hs_config = HypedsearchRunConfig.from_json(path=config.hs_config)
# logger.info("Getting missing outputs...")
missing_hybrid_outputs = hs_config.missing_hybrid_target_outputs
missing_native_outputs = hs_config.missing_native_target_outputs
# logger.info("Finished getting missing outputs")
_ = hs_config.print_missing_scan_info()

rule all:
    input:
        missing_hybrid_outputs + missing_native_outputs

rule run_hypedsearch:
    input: 
        mzml = lambda wildcards: next(
            Path(mzml) for mzml in hs_config.mzml_to_scans.keys() if Mzml.get_mzml_name(mzml) == wildcards.sample
        ),
    output:
        native_target = f"{hs_config.native_run_scan_results_dir}/{{sample}}.comet.{{scan}}-{{scan}}.target.txt",
        hybrid_target = f"{hs_config.hybrid_run_scan_results_dir}/{{sample}}.comet.{{scan}}-{{scan}}.target.txt"
    resources: 
        runtime = "10m"
    # benchmark:
    #     "logs/hypedsearch/{sample}.{scan}.log"
    singularity: 
        "docker://airikjohnson/hypedsearch:latest"
    shell:
        """
        python -m src.hypedsearch run-hypedsearch \
            --mzml {input.mzml} \
            --scan {wildcards.scan} \
            --config {config.hs_config}
        """