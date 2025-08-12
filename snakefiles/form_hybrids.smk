# NOTES
# Expected items in config file:
# 1. mzml_dir: directory containing mzML files
# 2. samples: list of MZML files from mzml_dir to include in analysis
# 3. hybrids_dir: directory to store hybrid JSON files as 
#     <hybrids_dir>/<sample>/<scan>.json
# 4. database: path to the k-mer database file from which we'll form hybrids from
# 5. fasta: path to the FASTA file to use for forming hybrids
# 6. precursor_mz_ppm_tol: PPM tolerance for precursor m/z
# 7. peak_to_ion_ppm_tol: PPM tolerance for peak to ion m/z$
# 8. job_name: name of the job to use in Snakemake

import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from types import SimpleNamespace
from src.mass_spectra import Mzml

##### Set other directory paths relative to the parent results directory
# Run on fewer scans if teesting
testing = bool(config.get("test", False))
num_scans_for_testing = 20
if testing:
    print(f"Running in testing mode on {num_scans_for_testing} scans.")

# Use python's SimpleNamespace to convert the config dict
# object into something that lets us refer to keys with dot notation
config = SimpleNamespace(**config)

# Get the path to all the expected hybrid JSONs
hybrid_outputs = []
for sample in config.samples:
    mzml = Mzml(path=f"{config.mzml_dir}/{sample}.mzML")
    sample_scan_paths = [f"{config.hybrids_dir}/{sample}/{scan}.json" for scan in mzml.scans]
    hybrid_outputs.extend(sample_scan_paths)

if testing:
    hybrid_outputs = hybrid_outputs[:num_scans_for_testing]

rule all:
    input: 
        hybrid_outputs

rule form_hybrids_for_scan:
    input: 
        mzml = f"{config.mzml_dir}/{{sample}}.mzML",
        db = config.database
    output:
        f"{config.hybrids_dir}/{{sample}}/{{scan}}.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        scan=lambda wildcards: int(wildcards.scan),
    # benchmark:
    #     "logs/snakemake/form_hybrids/test/{sample}_{scan}.txt"
    shell:
        """ 
        python -m src.hybrids_via_clusters \
        --mzml {input.mzml} \
        --scan {wildcards.scan} \
        --fasta {config.fasta} \
        --database {input.db} \
        --precursor_mz_ppm_tol {config.precursor_mz_ppm_tol} \
        --peak_to_ion_ppm_tol {config.peak_to_ion_ppm_tol} \
        --out_dir {config.hybrids_dir}/{wildcards.sample}
        """
        # # Use cProfile to profile the form_hybrids command
        # """
        # python -m cProfile -o logs/python/form_hybrids/test/{wildcards.sample}_{wildcards.scan}.prof \
        #     -m src.form_hybrids form-hybrids \
        #     --mzml {input.mzml} \
        #     --scan {wildcards.scan} \
        #     --fasta {config.fasta} \
        #     --db_path {input.db} \
        #     --precursor_mz_ppm_tol {config.precursor_mz_ppm_tol} \
        #     --peak_to_ion_ppm_tol {config.peak_to_ion_ppm_tol} \
        #     --out_dir {config.hybrids_dir}/{wildcards.sample}
        # """