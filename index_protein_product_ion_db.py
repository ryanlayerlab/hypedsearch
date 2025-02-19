from src.lookups.protein_product_ion_db import load_existing_protein_product_ion_db
import os
from src.erik_utils import setup_logger


logger = setup_logger()

db_path = "dbs/Uniprot_mouse.fasta_max_k=30_charges=(1, 2, 3).db"
size_gb = os.path.getsize(db_path) / (1024 ** 3)
logger.info(f"DB size BEFORE indexing = {size_gb}")

db = load_existing_protein_product_ion_db(db_path=db_path)
db.create_index_on_product_ion_mass()

size_gb = os.path.getsize(db_path) / (1024 ** 3)
logger.info(f"DB size AFTER indexing = {size_gb}")