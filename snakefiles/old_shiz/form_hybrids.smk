import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from src.mass_spectra import Mzml

# Set other directory paths relative to the parent results directory
config["native_run_dir"] = f"{config['parent_results_dir']}/native_run"
config["db_dir"] = f"{config['parent_results_dir']}/db"
config["hybrids_dir"] = f"{config['parent_results_dir']}/hybrids"
config["fasta_dir"] = f"{config['parent_results_dir']}/fastas"
config["hybrid_run_dir"] = f"{config['parent_results_dir']}/hybrid_run"


rule form_hybrids_for_scan:
    input: 
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML",
        db = f"{config['db_dir']}/top_{config['top_n_proteins']}_proteins.db"
    output:
        "{out_dir}/{sample}/{scan}.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        scan=lambda wildcards: int(wildcards.scan),
    shell:
        "python -m src.form_hybrids form-hybrids "
        "--mzml {input.mzml} "
        "--scan {wildcards.scan} "
        "--fasta {config[fasta]} "
        "--db_path {input.db} "
        "--precursor_mz_ppm_tol {config[precursor_mz_ppm_tol]} "
        "--peak_to_ion_ppm_tol {config[peak_to_ion_ppm_tol]} "
        "--out_dir {wildcards.out_dir}/{wildcards.sample} "