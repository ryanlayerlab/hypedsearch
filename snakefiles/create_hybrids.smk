configfile: "snakefiles/mouse_samples.yaml"
import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from src.mass_spectra import Mzml
from src.comet_utils import CometTxt
import numpy as np


mzml = Mzml(path="tests/data/mouse_spectra.mzML")
fasta = Path("tests/data/mouse_proteome.fasta").absolute()
db_path = Path("tests/data/mouse_top_10_prots.db").absolute()
config["mzml"] = str(mzml.path)
config["mzml_stem"] = mzml.path.stem
config["fasta"] = str(fasta)
config["db_path"] = str(db_path)
config["out_dir"] = f"tmp/{mzml.path.stem}"
config["num_psms"] = 10
# scans = mzml.scans

native_comet_run = Path("data/comet_run_1/BMEM_AspN_Fxn4/BMEM_AspN_Fxn4.txt").absolute()
psms = CometTxt(path=native_comet_run).read_psms()
scans = list(np.unique([psm.scan for psm in psms]))
scans = scans[:100]
# scans = list(range(1, 101))
# scans = scans[:10]
# scans = [7]

# outputs = [f"{config['out_dir']}/k={k}.json" for k in ks]

rule run_form_hybrids:
    input:
        [
            f"{config['out_dir']}/hybrids_{scan}.fasta"
            for scan in scans
        ] + 
        [
            f"{config['out_dir']}/hybrids_{scan}.json"
            for scan in scans
        ]
        
rule form_hybrids:
    output:
        hybrids_fasta="{out_dir}/hybrids_{scan}.fasta",
        hybrids_json="{out_dir}/hybrids_{scan}.json"
    shell:
        "python -m src.form_hybrids "
        "--mzml {config[mzml]} "
        "--scan {wildcards.scan} "
        "--db_path {config[db_path]} "
        "--fasta {config[fasta]} "
        "--out_dir {config[out_dir]} "
        "--include_fasta"


rule run_comet_on_hybrids:
    input:
        [
            f"{config['out_dir']}/comet.{scan}-{scan}.decoy.txt"
            for scan in scans
        ] + 
        [
            f"{config['out_dir']}/comet.{scan}-{scan}.target.txt"
            for scan in scans
        ]

rule hybrids_comet_run:
    input:
        hybrids_fasta="{out_dir}/hybrids_{scan}.fasta"
    output:
        "{out_dir}/comet.{scan}-{scan}.decoy.txt",
        "{out_dir}/comet.{scan}-{scan}.target.txt"
    shell:
        "/Users/erjo3868/repos/crux-4.3.Darwin.x86_64/bin/crux "
        "comet --parameter-file ./comet/crux.comet.params "
        "--decoy_search 2 "
        "--scan_range '{wildcards.scan} {wildcards.scan}' "
        "--output-dir {config[out_dir]} "
        "--num_output_lines {config[num_psms]} "
        # "--fileroot hybrids_{config[mzml_stem]} "
        "--overwrite T "
        "{config[mzml]} {input.hybrids_fasta}"

combined_targets = f"{config['out_dir']}/combined_targets.txt"
combined_decoys = f"{config['out_dir']}/combined_decoys.txt"

rule run_combine_targets_and_decoys:
    input:
        combined_targets,
        combined_decoys
    shell:
        "echo 'You now have: {input}'"

rule combine_targets_and_decoys:
    output: 
        combined_targets,
        combined_decoys,
    shell:
        """
        head -n 1 {config[out_dir]}/comet.7-7.target.txt > {output[0]}; 
        tail -n +2 -q {config[out_dir]}/comet.*-*.target.txt >> {output[0]}; 
        head -n 1 {config[out_dir]}/comet.7-7.decoy.txt > {output[1]}; 
        tail -n +2 -q {config[out_dir]}/comet.*-*.decoy.txt >> {output[1]}; 
        """
        
# rule run_comet_natively:
#     input:
#         hybrids_fasta=rules.form_hybrids.output.hybrids_fasta
#     output:
#         out_txt="{out_dir}/hybrids_{mzml_stem}.comet.{scan}-{scan}.txt"
#     shell:
#         "/Users/erjo3868/repos/crux-4.3.Darwin.x86_64/bin/crux "
#         "comet --parameter-file ./comet/crux.comet.params "
#         "--decoy_search 0 "
#         "--scan_range '{wildcards.scan} {wildcards.scan}' "
#         "--output-dir {config[out_dir]} "
#         "--num_output_lines {config[num_psms]} "
#         "--fileroot hybrids_{config[mzml_stem]} "
#         "--overwrite T "
#         "{config[mzml]} {input.hybrids_fasta}"



# rule all:
#     input:
#         outputs + [f"{config['out_dir']}/combined_ks.json"]

# rule create_kmer_to_prot_freq_pklz:
#     output:
#         "{out_dir}/k={k}.pklz"
#     shell:
#         "python -m src.peptides_and_ions get-kmer-info "
#         "-f {config[fasta]} "
#         "-k {wildcards.k} -o {config[out_dir]}/k={wildcards.k}.pklz"

# rule create_json:
#     input:
#         "{out_dir}/k={k}.pklz"
#     output:
#         "{out_dir}/k={k}.json"
#     script:
#         "analyze_kmer_fasta_data.py"

# rule combine_jsons:
#     input:
#         expand("{out_dir}/k={k}.json", out_dir=config["out_dir"], k=ks)
#     output:
#         "{out_dir}/combined_ks.json"
#     script:
#         "combine_kmer_jsons.py"