# Expected items in config file:
# 1. mzml_to_scans: dictionary mapping paths to mzML files to scan numbers
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

from types import SimpleNamespace
from src.mass_spectra import Mzml
from src.hypedsearch_utils import HypedsearchOutput
from dataclasses import asdict

# Use python's SimpleNamespace to convert the config dict
# object into something that lets us refer to keys with dot notation
config = SimpleNamespace(**config)

# Get expected outputs
expected_outputs = []
for mzml, scans in config.mzml_to_scans.items():
    if isinstance(scans, str):
        assert scans == "all", "Only 'all' scans are supported"
        scans = Mzml(path=mzml).scans

    for scan in scans:
        hs_outputs = HypedsearchOutput.get_hypedsearch_outputs(
            mzml=Path(mzml),
            scan=scan,
            out_dir=Path(config.out_dir),
        )
        outputs = [str(v) for v in asdict(hs_outputs).values()]
        expected_outputs.extend(outputs)

rule all:
    input: 
        expected_outputs

rule run_hypedsearch_on_spectrum:
    input:
        mzml=lambda wildcards: next(
            Path(m) for m in config.mzml_to_scans.keys() if Path(m).stem == wildcards.sample
        ),
        fasta = config.fasta,
        crux_comet_params = config.crux_comet_params,
        database = config.database,
    output:
        native_target = "{out_dir}/native_{sample}.comet.{scan}-{scan}.target.txt",
        native_decoy = "{out_dir}/native_{sample}.comet.{scan}-{scan}.decoy.txt",
        hybrid_target = "{out_dir}/hybrid_{sample}.comet.{scan}-{scan}.target.txt",
        hybrid_decoy = "{out_dir}/hybrid_{sample}.comet.{scan}-{scan}.decoy.txt",
        hybrids = "{out_dir}/{sample}_scan={scan}.json",
    shell:
        """
        python -m src.hypedsearch_utils \
            --mzml {input.mzml} \
            --scan {wildcards.scan} \
            --database {input.database} \
            --fasta {input.fasta} \
            --crux_comet_params {input.crux_comet_params} \
            --out_dir {wildcards.out_dir}
        """

