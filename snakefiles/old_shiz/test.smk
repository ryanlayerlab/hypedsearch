configfile: "snakefiles/mouse_samples.yaml"

from pathlib import Path

rule print_config:
    run:
        print("Config dictionary:", config)
        # print("MZMLs:", path_mzmls)

rule all:
    input:
        expand(
            "{out_dir}/native/{sample}.comet.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        ),
        expand(
            "{out_dir}/native/{sample}.scans.txt", 
            out_dir=config["out_dir"], 
            sample=config["samples"]
        )



rule run_comet_natively:
    input: 
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
    output:
        comet_out = f"{config['out_dir']}/native/{{sample}}.comet.txt"
    shell:
        "/Users/erjo3868/repos/crux-4.3.Darwin.x86_64/bin/crux "
        "comet --parameter-file ./comet/crux.comet.params "
        "--decoy_search 0 "
        "--output-dir {config[out_dir]}/native "
        "--num_output_lines {config[num_psms]} "
        "--fileroot {wildcards.sample} "
        "--overwrite T "
        "{input.mzml} {config[fasta]} && "
        "rm {config[out_dir]}/native/{wildcards.sample}.comet.*.txt"


rule cut_scans_with_results:
    input:
        comet_out = f"{config['out_dir']}/native/{{sample}}.comet.txt"
    output:
        scans_txt = f"{config['out_dir']}/native/{{sample}}.scans.txt"
    shell:
        "tail -n +2 {input.comet_out} | cut -f 1 | uniq > {output.scans_txt}"

rule form_hybrids:
    input: 
        scans_txt = f"{config['out_dir']}/native/{{sample}}.scans.txt",
        mzml = f"{config['mzml_dir']}/{{sample}}.mzML"
    output:
        hybrids_fasta = f"{config['out_dir']}/hybrids/{{sample}}/.scans.txt""{out_dir}/hybrids_{scan}.fasta",
        hybrids_json="{out_dir}/hybrids_{scan}.json"
    shell:
        "python -m src.form_hybrids "
        "--mzml {input.mzml} "
        "--scans_path {input.scans_txt} "
        "--db_path {config[db_path]} "
        "--fasta {config[fasta]} "
        "--out_dir {config[out_dir]}/hybrids/{wildcards.sample} "
        "--include_fasta"

rule run_comet_on_hybrids:
    input:
        