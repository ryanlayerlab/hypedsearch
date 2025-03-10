import sys
from pathlib import Path
from time import time
from collections import defaultdict

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from src.utils import decompress_and_unpickle, pickle_and_compress
from src.utils import log_params, pickle_and_compress, setup_logger, get_time_in_diff_units
from src.constants import RESULTS_DIR

logger = setup_logger()

testing = False

@log_params
def main(uniq_peptides_path):
    logger.info(f"Loading uniq_peptides dict...")
    load_start = time()
    uniq_peptides = decompress_and_unpickle(uniq_peptides_path)
    total_time = get_time_in_diff_units(time()-load_start)
    logger.info(f"Loading uniq_peptides dict took {total_time}")
    logger.info(f"There are {len(uniq_peptides)} unique peptides")
    
    logger.info("Getting dictionary of peptide length vs num proteins it appears in")
    num_proteins_by_len = defaultdict(list)
    for seq, proteins in uniq_peptides.items():
        num_proteins_by_len[len(seq)].append(len(proteins))

    # Save the data
    file_name = f"{uniq_peptides_path.stem}_num_proteins_by_len.pkl"
    if testing:
        file_name = f"testing/{file_name}"
    out_path = RESULTS_DIR / file_name
    logger.info(f"Saving data to {out_path}")

    # Option 1 - compressed pickle
    pickle_and_compress(obj=num_proteins_by_len, file_path=out_path)

    # # Option 2 - joblib
    # dump(uniq_peptides, out_path, compress=3)

    logger.info("Finished saving the data. Ending!")


if __name__ == "__main__":
    if testing:
        uniq_peptides_path = Path("results/testing/test_uniqPeptides_minK=1_maxK=30.pkl")

    else:
        uniq_peptides_path = Path("results/Uniprot_mouse_uniqPeptides_minK=1_maxK=30.pkl")
    main(uniq_peptides_path)
