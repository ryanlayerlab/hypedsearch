import sys
from pathlib import Path
from types import SimpleNamespace

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from src.mass_spectra import Mzml

# Set other directory paths relative to the parent results directory
config["native_run_dir"] = f"{config['parent_results_dir']}/native_run"
config["db_dir"] = f"{config['parent_results_dir']}/db"
config["hybrids_dir"] = f"{config['parent_results_dir']}/hybrids"
config["fasta_dir"] = f"{config['parent_results_dir']}/fastas"
config["hybrid_run_dir"] = f"{config['parent_results_dir']}/hybrid_run"
config["native_hybrid_run_dir"] = f"{config['parent_results_dir']}/native_hybrid_run"

# Make directories if they don't exist
for dir_name in ["native_run_dir", "db_dir", "hybrids_dir", "fasta_dir", "hybrid_run_dir", "native_hybrid_run_dir"]:
    Path(config[dir_name]).mkdir(parents=True, exist_ok=True)

# Get all scans for each sample
hybrid_outputs = []
hybrid_outputs_by_sample = {}
for sample in config["samples"]:
    mzml = Mzml(path=f"{config['mzml_dir']}/{sample}.mzML")
    sample_scan_paths = [f"{config['hybrids_dir']}/{sample}/{scan}.json" for scan in mzml.scans]
    hybrid_outputs.extend(sample_scan_paths)
    hybrid_outputs_by_sample[sample] = sample_scan_paths

# Use python's SimpleNamespace to convert the config dict
# object into something that let's me refer to keys with dot notation
config = SimpleNamespace(**config)

rule all:
    input:
        # Native Comet + assign-confidence run
        expand(
            "{out_dir}/{sample}.comet.target.txt",
            out_dir=config.native_run_dir,
            sample=config.samples
        ),
        expand(
            "{out_dir}/{sample}.comet.decoy.txt",
            out_dir=config.native_run_dir,
            sample=config.samples
        ),
        expand(
            "{out_dir}/assign-confidence.target.txt",
            out_dir=config.native_run_dir,
        ),
        # Protein and prefix abundances
        f"{config.db_dir}/top_{config.top_n_proteins}_proteins.txt",
        f"{config.db_dir}/protein_abundances.png",
        # f"{config.native_run_dir}/prefix_abundances.json",
        # Create database from (seemingly) most abundant proteins
        f"{config.db_dir}/top_{config.top_n_proteins}_proteins.db",
        # Form hybrids
        
        # expand(
        #     "{out_dir}/{sample}.json",
        #     out_dir=config["hybrids_dir"],
        #     sample=config["samples"]
        # )
        # Method 1: form hybrids for all scans in each mzML
        # expand(
        #     "{out_dir}/{sample}/done.txt",
        #     out_dir=config["hybrids_dir"],
        #     sample=config["samples"]
        # )
        # Method 2: form hybrids for each scan one at a time -- more easily parallelizable
        # at the cost of having to load the MZML and locate the scan each time
        hybrid_outputs,

        # Create FASTA file for each sample
        # expand(
        #     "{out_dir}/{sample}.fasta",
        #     out_dir=config["fasta_dir"],
        #     sample=config["samples"]
        # )

        # Comet run with hybrids
        expand(
            "{out_dir}/{sample}.comet.target.txt",
            out_dir=config.hybrid_run_dir,
            sample=config.samples
        ),
        expand(
            "{out_dir}/{sample}.comet.decoy.txt",
            out_dir=config.hybrid_run_dir,
            sample=config.samples
        ),
        expand(
            "{out_dir}/assign-confidence.target.txt",
            out_dir=config.hybrid_run_dir,
        ),

        # # Comet run with native+hybrids
        # expand(
        #     "{out_dir}/{sample}.comet.target.txt",
        #     out_dir=config["native_hybrid_run_dir"],
        #     sample=config["samples"]
        # ),
        # expand(
        #     "{out_dir}/{sample}.comet.decoy.txt",
        #     out_dir=config["native_hybrid_run_dir"],
        #     sample=config["samples"]
        # ),
        # expand(
        #     "{out_dir}/assign-confidence.target.txt",
        #     out_dir=config["native_hybrid_run_dir"],
        # ),



rule comet_separate_decoy_search:
    input:
        mzml = f"{config.mzml_dir}/{{sample}}.mzML",
        fasta = lambda wildcards: 
            config.fasta if wildcards.out_dir == config.native_run_dir else config.fasta_dir + f"/{wildcards.sample}.fasta"
    output:
        target_out = "{out_dir}/{sample}.comet.target.txt",
        decoy_out = "{out_dir}/{sample}.comet.decoy.txt"
    shell:
        "./comet/crux-4.3.Darwin.x86_64/bin/crux comet "
        "--parameter-file {config.crux_comet_params} "
        "--decoy_search 2 "
        "--output-dir {wildcards.out_dir} "
        "--fileroot {wildcards.sample} "
        "--verbosity 60 "
        "--overwrite T "
        "--spectrum_batch_size 500 "
        "{input.mzml} {input.fasta}; "
        "rm {wildcards.out_dir}/{wildcards.sample}.comet.params.txt; " # remove params
        "rm {wildcards.out_dir}/{wildcards.sample}.comet.log.txt" # remove log file

rule assign_confidence:
    input:
        lambda wildcards: [
            f"{wildcards.out_dir}/{sample}.comet.target.txt"
            for sample in config.samples 
        ]
    output:
        "{out_dir}/assign-confidence.target.txt"
    shell:
        "./comet/crux-4.3.Darwin.x86_64/bin/crux assign-confidence "
        "--output-dir {wildcards.out_dir} "
        "--overwrite T "
        "--list-of-files T "
        "{input}; "
        "rm {wildcards.out_dir}/assign-confidence.params.txt; " # remove params
        "rm {wildcards.out_dir}/assign-confidence.log.txt" # remove log file"


rule get_most_common_proteins:
    input:
        f"{config.native_run_dir}/assign-confidence.target.txt"
    output:
        f"{config.db_dir}/top_{config.top_n_proteins}_proteins.txt",
        f"{config.db_dir}/protein_abundances.png"
    shell:
        "python -m src.protein_abundance protein-abundances "
        "--comet_results_dir {config.native_run_dir} "
        "--q_value_threshold {config.q_value_threshold} "
        "--top_n_proteins {config.top_n_proteins} "
        "--out_path {config.db_dir}/top_{config.top_n_proteins}_proteins.txt " 

rule get_prefix_abundances:
    input: 
        f"{config.native_run_dir}/assign-confidence.target.txt"
    output:
        f"{config.native_run_dir}/prefix_abundances.json"
    shell:
        "python -m src.protein_abundance prefix-abundances "
        "--comet_results_dir {config.native_run_dir} "
        "--q_value_threshold {config.q_value_threshold} " 
        "--out_path {config.native_run_dir}/prefix_abundances.json"

rule create_db:
    input: 
        top_proteins_txt = f"{config.db_dir}/top_{config.top_n_proteins}_proteins.txt"
    output:
        top_proteins_db = f"{config.db_dir}/top_{config.top_n_proteins}_proteins.db"
    shell:
        "python -m src.create_db "
        "--db_path {output.top_proteins_db} "
        "--proteins {input.top_proteins_txt} "
        "--fasta_path {config.fasta}"

rule form_hybrids_for_scan:
    input: 
        mzml = f"{config.mzml_dir}/{{sample}}.mzML",
        db = f"{config.db_dir}/top_{config.top_n_proteins}_proteins.db"
    output:
        "{out_dir}/{sample}/{scan}.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        scan=lambda wildcards: int(wildcards.scan),
    shell:
        "python -m src.form_hybrids form-hybrids "
        "--mzml {input.mzml} "
        "--scan {wildcards.scan} "
        "--fasta {config.fasta} "
        "--db_path {input.db} "
        "--precursor_mz_ppm_tol {config.precursor_mz_ppm_tol} "
        "--peak_to_ion_ppm_tol {config.peak_to_ion_ppm_tol} "
        "--out_dir {wildcards.out_dir}/{wildcards.sample} "

rule create_hybrids_fasta:
    input:
        hybrid_outputs = lambda wc: hybrid_outputs_by_sample[wc.sample]
    output:
        # The FASTA file is temporary because it can be large (>1 GB) and it contains
        # the native FASTA sequences and the hybrid sequences both of which are elsewhere
        # so there's no need to keep it around
        temp(f"{config.fasta_dir}/{{sample}}.fasta")
    shell:
        "python -m src.form_hybrids hybrids-to-fasta "
        "--hybrids_path {config.hybrids_dir}/{wildcards.sample} "
        "--old_fasta {config.fasta} " # include native FASTA for decoy database search purposes
        "--new_fasta_path {config.fasta_dir}/{wildcards.sample}.fasta"

