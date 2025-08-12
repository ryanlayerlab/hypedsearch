import sys
from pathlib import Path

repo_dir = Path(workflow.basedir).absolute().parent
sys.path.append(str(repo_dir))

from types import SimpleNamespace

fasta_stem = Path(config["fasta"]).stem
# Use python's SimpleNamespace to convert the config dict
# object into something that let's me refer to keys with dot notation
config = SimpleNamespace(**config)

kmer_to_protein_json = f"{config.out_dir}/{fasta_stem}.pklz"
kmer_db = f"{config.out_dir}/{fasta_stem}.db"

rule all:
    input:
        kmer_to_protein_json,
        kmer_db


rule create_kmer_to_protein_map:
    input:
        fasta = config.fasta
    output:
        kmer_to_protein_json = kmer_to_protein_json,
    # shell:
    #     "python -m src.peptides_and_ions create-kmer-to-protein-map "
    #     "--fasta {config.fasta} "
    #     "--protein_attr name "
    #     "--output_path {config.out_dir}/{fasta_stem}.json"
    run:
        from src.utils import setup_logger
        from src.create_db import KmerToProteinMap
        setup_logger()
        kmer_to_protein_map = KmerToProteinMap.create(
            fasta=input.fasta,
            protein_attr="name",
        )
        kmer_to_protein_map.save(output.kmer_to_protein_json)



rule create_db:
    input:
        kmer_to_protein_json = kmer_to_protein_json,
    output:
        kmer_db = kmer_db
    run:
        from src.utils import setup_logger
        from src.create_db import KmerDatabase
        setup_logger()
        kmer_db = KmerDatabase.create_db(
            db_path=output.kmer_db,
            kmer_to_protein_map=input.kmer_to_protein_json,
        )



