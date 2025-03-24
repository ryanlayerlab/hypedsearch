import sys
from pathlib import Path
from time import time

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from src.constants import MOUSE_PROTEOME, RESULTS_DIR
from src.peptides_and_ions import Peptide, get_uniq_kmer_to_protein_map, write_fasta
from src.utils import log_params, pickle_and_compress, setup_logger

logger = setup_logger()

testing = False


@log_params
def main(fasta_path: str, min_k: int, max_k: int):
    proteins = Peptide.from_fasta(fasta_path=fasta_path)
    uniq_peptides = get_uniq_kmer_to_protein_map(min_k=1, max_k=30, proteins=proteins)
    file_name = f"{fasta_path.stem}_uniqPeptides_minK={min_k}_maxK={max_k}.pkl"
    out_path = RESULTS_DIR / file_name
    pickle_and_compress(obj=uniq_peptides, file_path=out_path)


if __name__ == "__main__":
    if testing:
        fasta_path = Path("tmp/test.fasta")
        proteins = Peptide.from_fasta(fasta_path=MOUSE_PROTEOME)
        write_fasta(peptides=proteins[:10], output_path=fasta_path)
    else:
        fasta_path = MOUSE_PROTEOME
    main(fasta_path=fasta_path, min_k=1, max_k=30)
