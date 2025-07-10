import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

# Set other directory paths relative to the parent results directory
config["native_run_dir"] = f"{config['parent_results_dir']}/native_run"

rule all:
    input:
        # Native Comet + assign-confidence run
        expand(
            "{out_dir}/{sample}.comet.target.txt",
            out_dir=config["native_run_dir"],
            sample=config["samples"]
        ),
        expand(
            "{out_dir}/{sample}.comet.decoy.txt",
            out_dir=config["native_run_dir"],
            sample=config["samples"]
        ),
        # expand(
        #     "{out_dir}/assign-confidence.target.txt",
        #     out_dir=config["native_run_dir"],
        # )

rule comet_via_singularity:
    input:
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML",
        fasta = lambda wildcards: 
            config["fasta"] if wildcards.out_dir == config["native_run_dir"] else config["fasta_dir"] + f"/{wildcards.sample}.fasta"
    output:
        target_out = "{out_dir}/{sample}.comet.target.txt",
        decoy_out = "{out_dir}/{sample}.comet.decoy.txt"
    params:
        decoy_search = 2
    shell:
        "singularity exec "
        "--bind {input.mzml}:/data/mzml.mzML "
        "--bind {input.fasta}:/data/fasta.fasta "
        "--bind {config.crux_comet_params}:/data/crux.comet.params "
        "--bind {wildcards.out_dir}:/results "
        "{config.singularity_image} "
        "comet "
        "--parameter-file /data/crux.comet.params "
        "--decoy_search {params.decoy_search} "
        "--output-dir /results "
        "--fileroot {wildcards.sample} "
        "--verbosity 60 "
        "--overwrite T "
        "--spectrum_batch_size 500 "
        " {input.fasta}; "
        "rm /results/{wildcards.sample}.comet.params.txt; " # remove params
        "rm /results/{wildcards.sample}.comet.log.txt" # remove log file

rule assign_confidence:
    input:
        lambda wc: [
            f"{wc.out_dir}/{sample}.comet.target.txt"
            for sample in config["samples"] 
        ]
    output:
        "{out_dir}/assign-confidence.target.txt"
    shell:
        "singularity exec {config.singularity_image} "
        "crux assign-confidence "
        "--output-dir {wildcards.out_dir} "
        "--overwrite T "
        "--list-of-files T "
        "{input}; "
        "rm {wildcards.out_dir}/assign-confidence.params.txt; " # remove params
        "rm {wildcards.out_dir}/assign-confidence.log.txt" # remove log file"