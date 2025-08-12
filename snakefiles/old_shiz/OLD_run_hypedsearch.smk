
config["out_dir"] = Path(config["out_dir"]).absolute()
config["mzml_dir"] = Path(config["mzml_dir"]).absolute()
config["db_path"] = Path(config["db_path"]).absolute()

import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))


rule run_all:
    input:
        # Form hybrids
        expand(
            "{out_dir}/hybrids/{sample}.json", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        ),
        # Native Comet run
        expand(
            "{out_dir}/native_run/{sample}.comet.target.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        ),
        expand(
            "{out_dir}/native_run/{sample}.comet.decoy.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )



# rule run_hypedsearch_on_mzml:
#     input:
#         expand(
#             "{out_dir}/{sample}/done.txt", 
#             out_dir=config["out_dir"], 
#             sample=config["samples"]
#         )

# rule hypedsearch_on_mzml:
#     input:
#         mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
#     output:
#         f"{config['out_dir']}/{{sample}}/done.txt"
#     shell:
#         "python -m src.run_hypedsearch -m {input.mzml} "
#         "-d {config[db_path]} -f {config[fasta]} -o {config[out_dir]} "
#         "-n {config[num_psms]} "
#         # "-t" # for testing

rule run_form_hybrids_for_mzmls:
    input:
        expand(
            "{out_dir}/hybrids/{sample}.json", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )

rule form_hybrids_for_mzml:
    input:
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
    output:
        json = f"{config['out_dir']}/hybrids/{{sample}}.json"
    shell:
        "python -m src.form_hybrids "
        "-m {input.mzml} "
        "-f {config[fasta]} "
        "-d {config[db_path]} "
        "-pmpt {config[precursor_mz_ppm_tol]} "
        "-pipt {config[peak_to_ion_ppm_tol]} "
        "-o {output.json}"

rule run_postprocess_hybrids:
    input:
        expand(
            "{out_dir}/hybrids/{sample}.json", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )

    output:
        expand(
            "{out_dir}/hybrids/all_hybrids.json", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )


rule run_native_comet:
    input:
        expand(
            "{out_dir}/native_run/{sample}.comet.target.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        ),
        expand(
            "{out_dir}/native_run/{sample}.comet.decoy.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )


rule comet:
    input:
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
    output:
        target_out = f"{config['out_dir']}/{{sample}}.comet.target.txt",
        decoy_out = f"{config['out_dir']}/{{sample}}.comet.decoy.txt"
    
rule native_comet:
    input:
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
    output:
        target_out = f"{config['out_dir']}/native_run/{{sample}}.comet.target.txt",
        decoy_out = f"{config['out_dir']}/native_run/{{sample}}.comet.decoy.txt"
    shell:
        "/Users/erjo3868/repos/crux-4.3.Darwin.x86_64/bin/crux comet "
        "--parameter-file {config[comet_params]} "
        "--decoy_search 2 "
        "--output-dir {config[out_dir]}/native_run "
        "--fileroot {wildcards.sample} "
        "--verbosity 60 "
        # "--overwrite T "
        "--scan_range '0 100' " # for testing
        "--spectrum_batch_size 500 "
        "{input.mzml} {config[fasta]}; "
        "rm {config[out_dir]}/native_run/{wildcards.sample}.comet.params.txt; "
        "rm {config[out_dir]}/native_run/{wildcards.sample}.comet.log.txt"

# rule run_hybrid_comet:
#     input:
#         expand(
#             "{out_dir}/hybrid_run/{sample}.comet.target.txt", 
#             out_dir=config["out_dir"], 
#             sample=config["samples"]
#         )

# rule hybrid_comet:
#     input:
#         mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
#     output:
#         target_out = f"{config['hybrid_run_dir']}/{{sample}}.comet.target.txt",
#         decoy_out = f"{config['hybrid_run_dir']}/{{sample}}.comet.decoy.txt"
#     shell:
#         "/Users/erjo3868/repos/crux-4.3.Darwin.x86_64/bin/crux comet "
#         "--parameter-file {config[comet_params]} "
#         "--decoy_search 2 "
#         "--output-dir {config[hybrid_run_dir]} "
#         "--fileroot {wildcards.sample} "
#         "--verbosity 60 "
#         "--spectrum_batch_size 500 "
#         "--overwrite T "
#         "{input.mzml} {config[fasta_dir]}/{wildcards.sample}_hybrids.fasta"