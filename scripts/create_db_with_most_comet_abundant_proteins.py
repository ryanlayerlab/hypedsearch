import sys
from pathlib import Path

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from time import time

from src.constants import MOUSE_PROTEOME, RESULTS_DIR, IonTypes
from src.peptides_and_ions import get_proteins_by_name
from src.protein_product_ion_database import (
    create_and_populate_protein_and_product_ion_database,
)
from src.utils import decompress_and_unpickle, pickle_and_compress, setup_logger

logger = setup_logger()

protein_abundances = decompress_and_unpickle("results/protein_abundances.pkl")

# charges, ion_types = [1, 2, 3], [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
top_n_most_abundant_proteins = 3
fasta_path = MOUSE_PROTEOME
out_path = (
    RESULTS_DIR / f"db_{top_n_most_abundant_proteins}_most_abundant_comet_prots.db"
)

protein_names = list(protein_abundances[:top_n_most_abundant_proteins]["full_name"])
# protein_names, len(protein_names)

db_proteins = get_proteins_by_name(fasta_path=fasta_path, protein_names=protein_names)
start_time = time()
create_and_populate_protein_and_product_ion_database(
    # charges=charges, ion_types=ion_types,
    peptides=db_proteins,
    db_path=out_path,
)
# pickle_and_compress(obj=db, file_path=out_path)
logger.info(f"Took {round(time() - start_time, 2)} seconds")
