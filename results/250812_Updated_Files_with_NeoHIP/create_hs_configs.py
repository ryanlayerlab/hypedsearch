import sys
from pathlib import Path

from src.hypedsearch import HybridRunConfig
from src.mass_spectra import Mzml
from src.utils import setup_logger

logger = setup_logger()

mzml_dir = "data/spectra/human_samples/250812_Updated_Files_with_NeoHIP"

mzml_sets = {
    "3A": [
        f"{mzml_dir}/HuIslet_240121_3A_Baseline_AspN_Fxn{fxn}.mzML"
        for fxn in [3, 4, 5, 6, 7, 8]
    ],
    "3B": [
        f"{mzml_dir}/HuIslet_240121_3B_IL1beta_AspN_Fxn{fxn}.mzML"
        for fxn in [3, 4, 5, 6, 7, 8]
    ],
    "2A": [
        f"{mzml_dir}/HuIslets_250525_MD_2A_Baseline_NoTreatment_AspN_Fxn{fxn}.mzML"
        for fxn in [2, 3, 4, 5, 6, 7]
    ],
    "2B": [
        f"{mzml_dir}/HuIslets_250525_MD_2B_IL1B_24hr_AspN_Fxn{fxn}.mzML"
        for fxn in [2, 3, 4, 5, 6, 7]
    ],
    "1A": [
        f"{mzml_dir}/HuIslets_250620_1A_Baseline_20hrNoTreatment_AspN_Fxn{fxn}.mzML"
        for fxn in [5, 6, 7, 8]
    ],
    "1B": [
        f"{mzml_dir}/HuIslets_250620_1B_IL1B_20hr_AspN_Fxn{fxn}.mzML"
        for fxn in [4, 5, 6, 7, 8]
    ],
}

# The argument is the second item in sys.argv (first is the script name)
name = sys.argv[1]

parent_out_dir = "results/250812_Updated_Files_with_NeoHIP"
logger.info(f"Creating Hypedsearch config for set {name}...")
mzmls = mzml_sets[name]
mzml_to_scans = {}
for idx, mzml in enumerate(mzmls):
    logger.info(f"Loading MZML {idx+1} of {len(mzmls)}: {mzml}...")
    mzml_to_scans[mzml] = Mzml(mzml=mzml).scans
    # For testing
    # if idx >= 1:
    #     break

logger.info(f"Creating and saving config...")
results_dir = f"{parent_out_dir}/{name}/top35prots"
config = HybridRunConfig(
    mzml_to_scans=mzml_to_scans,
    out_dir=f"{results_dir}/scan_results",
    fasta="fastas/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fasta",
    kmer_db=f"{results_dir}/{name}_kmers.db",
    kmer_to_protein_map=f"{results_dir}/{name}_kmer_to_proteins.json",
    crux_path="/usr/local/bin/crux",
    peak_to_ion_ppm_tol=20,
    precursor_mz_ppm_tol=20,
    crux_comet_params=f"{parent_out_dir}/crux.comet.params",
)
node_data_dir = Path("/cache")
out_config_path = node_data_dir / f"{name}.config"
out_config_path = Path(sys.argv[2])
node_config = config.create_config_for_fiji_node(
    node_data_dir=node_data_dir, out_path=out_config_path
)
