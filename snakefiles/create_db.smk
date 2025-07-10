
rule all:
    input:
        f"{config['db_path']}/top_{config['top_n_proteins']}_proteins.txt",
        f"{config['db_path']}/protein_abundances.png"

rule get_most_common_proteins:
    input:
        "{out_dir}/assign-confidence.target.txt"
    output:
        f"{{db_dir}}/top_{config['top_n_proteins']}_proteins.txt",
        f"{config['db_dir']}/protein_abundances.png"
    shell:
        "python -m src.protein_abundance "
        "--comet_results_dir {config[out_dir]} "
        "--q_value_threshold {config[q_value_threshold]} "
        "--top_n_proteins {config[top_n_proteins]} "
        "--out_path {config[db_dir]}/top_{config[top_n_proteins]}_proteins.txt "  