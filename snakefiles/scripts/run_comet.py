import sys
from pathlib import Path
from types import SimpleNamespace

from snakemake.script import snakemake

from src.crux import CometConfig, Crux
from src.mass_spectra import Mzml
from src.utils import setup_logger

# config = SimpleNamespace(**snakemake.config)
mzml = Path(snakemake.input.mzml)
out_file = Path(snakemake.output.target)
scan = int(snakemake.wildcards.scan)
sample = snakemake.wildcards.sample
out_dir = Path(snakemake.wildcards.out_dir)

# Setup logger
# log_file = Path(
#     f"logs/{Mzml.get_mzml_name(snakemake.input.mzml)}_{snakemake.wildcards.scan}.log"
# )
logger = setup_logger()
logger.debug("Starting Comet search via run_comet.py script")
logger.debug(
    f"Running with params:\n\tmzml={mzml},\n\tscan={scan},\n\tsample={sample},\n\t"
    f"out_dir={out_dir},\n\t Snakemake config: {snakemake.config}"
)
crux = Crux()
crux.run_comet(
    mzml=snakemake.input.mzml,
    fasta=snakemake.config.fasta,
    crux_comet_params=snakemake.config.crux_comet_params,
    decoy_search=snakemake.config.decoy_search,
    out_dir=snakemake.config.out_dir,
    file_root=snakemake.wildcards.sample,
    scan_min=int(snakemake.wildcards.scan),
    scan_max=int(snakemake.wildcards.scan),
    num_threads=int(snakemake.config.num_threads),
)
