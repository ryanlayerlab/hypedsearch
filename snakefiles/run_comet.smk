
import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from types import SimpleNamespace
from src.mass_spectra import Mzml
import time


config["mzml_dir"] = "/localscratch/spectra"
config["samples"] = [   
  "HuIslets_A5_Sample2A_AspN_Fxn5_250520_MD",
  "HuIslets_B5_Sample2B_AspN_Fxn5_250520_MD",
]
config["hybrids_dir"] = "/localscratch/hybrids"
config["fasta"] = "/localscratch/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fasta"
config["crux_comet_params"] = "/localscratch/crux.comet.params"
config["fasta_dir"] = "/localscratch/fastas"
config["comet_results_dir"] = "/localscratch/comet"

# Use python's SimpleNamespace to convert the config dict
# object into something that lets us refer to keys with dot notation
config = SimpleNamespace(**config)

# Get the path to all the expected hybrid JSONs
print("Gathering outputs...")
t0 = time.perf_counter()
comet_outputs = []
for sample in config.samples:
    mzml = Mzml(path=f"{config.mzml_dir}/{sample}.mzML")
    comet_outputs.extend(
        [f"{config.comet_results_dir}/{sample}.comet.{scan}-{scan}.target.txt" for scan in mzml.scans]
    )
# comet_outputs = comet_outputs[:100]
# comet_outputs = [optional(p) for p in comet_outputs]
print(f"Finished gathering outputs. Took {round(time.perf_counter() - t0, 2)} seconds.")

rule all:
    input: 
        comet_outputs

rule create_hybrids_fasta:
    input:
        fasta = config.fasta,
        hybrid_json = f"{config.hybrids_dir}/{{sample}}_{{scan}}.json"
    output:
        hybrid_fasta = temp(
            f"{config.fasta_dir}/{{sample}}_{{scan}}.fasta"
        )
    shell:
        """
        python -m src.form_hybrids hybrids-to-fasta \
            --hybrids_path {input.hybrid_json} \
            --new_fasta_path {output.hybrid_fasta} \
            --old_fasta {input.fasta}
        """

rule comet_via_singularity_on_scan:
    input:
        mzml = f"{config.mzml_dir}/{{sample}}.mzML",
        fasta = config.fasta,
        hybrid_fasta = f"{config.fasta_dir}/{{sample}}_{{scan}}.fasta"
    output:
        target_out = "{out_dir}/{sample}.comet.{scan}-{scan}.target.txt",
        decoy_out = "{out_dir}/{sample}.comet.{scan}-{scan}.decoy.txt"
    # params:
        # decoy_search = 2
    shell:
        """
        tmpdir=$(mktemp -d)
        set +e  # turn off 'exit on error'
        singularity exec \
            --bind {input.mzml}:/data/mzml.mzML \
            --bind {input.hybrid_fasta}:/data/fasta.fasta \
            --bind {config.crux_comet_params}:/data/crux.comet.params \
            --bind $tmpdir:/results \
            hypedsearch_latest.sif \
            bash -c "crux comet \
            --parameter-file /data/crux.comet.params \
            --num_threads 1 \
            --decoy_search 2 \
            --fileroot {wildcards.sample} \
            --scan_range '{wildcards.scan} {wildcards.scan}' \
            --output-dir /results \
            --verbosity 60 \
            --overwrite T \
            --spectrum_batch_size 500 \
            /data/mzml.mzML /data/fasta.fasta"
        
        set -e  # restore default behavior if needed

        # Copy the expected outputs where they should go and create empty ones if Comet didn't
        # spit anything out
        target_file=$tmpdir/{wildcards.sample}.comet.{wildcards.scan}-{wildcards.scan}.target.txt
        if [ -f "$target_file" ]; then
            cp "$target_file" {wildcards.out_dir}
        else
            touch {output.target_out}
        fi

        decoy_file=$tmpdir/{wildcards.sample}.comet.{wildcards.scan}-{wildcards.scan}.decoy.txt
        if [ -f "$decoy_file" ]; then
            cp "$decoy_file" {wildcards.out_dir}
        else
            touch {output.decoy_out}
        fi
        """
        
        # "singularity exec "
        # "--bind {input.mzml}:/data/mzml.mzML "
        # "--bind {input.fasta}:/data/fasta.fasta "
        # "--bind {config.crux_comet_params}:/data/crux.comet.params "
        # "--bind {wildcards.out_dir}:/results "
        # "{config.singularity_image} "
        # "bash -c '"
        # "crux comet "
        # "--parameter-file /data/crux.comet.params "
        # "--decoy_search {params.decoy_search} "
        # "--output-dir /results "
        # "--fileroot {wildcards.sample} "
        # "--verbosity 60 "
        # "--overwrite T "
        # "--spectrum_batch_size 500 "
        # "/data/mzml.mzML /data/fasta.fasta "
        # "&& rm /results/{wildcards.sample}.comet.params.txt "
        # "&& rm /results/{wildcards.sample}.comet.log.txt'"