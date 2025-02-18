from pathlib import Path

# from src.erik_utils import file_exists

# Paths
GIT_REPO_DIR = Path(__file__).parents[1]
ROOT_DIR = Path(__file__).parents[2]
TEST_DIR = GIT_REPO_DIR / "tests"
LOGS_DIR = GIT_REPO_DIR / "logs"
PLOTS_DIR = GIT_REPO_DIR / "plots"
DATA_DIR = GIT_REPO_DIR / "data"
COMET_RUN_1_DIR = DATA_DIR / "comet_run_1"
COMET_RUN_2_DIR = DATA_DIR / "comet_run_2"
HS_DIR = DATA_DIR / "hs"
SPECTRA_DIR = DATA_DIR / "spectra"
COMET_EXECUTABLE = GIT_REPO_DIR / f"comet/comet.macos.exe"
COMET_PARAMS = GIT_REPO_DIR / f"comet/comet.params"

FASTAS_DIR = GIT_REPO_DIR / "fastas"
MOUSE_PROTEOME = FASTAS_DIR / "Uniprot_mouse.fasta"

# Strings
MASS = "mass"
KMERS = "kmers"
SEQ = "sequence"
COUNT = "count"
BOTH = "both"
KMER_AND_MASS = "kmer+mass"
KMER = "kmer"
SAMPLE = "sample"
SPECTRUM_ID = "spectrum_id"
SCAN = "scan"
PLAIN_PEPTIDE = "plain_peptide"
CHARGE = "charge"
INCLUSIVE_START = "inclusive_start"
I_START = "inclusive_start"
EXCLUSIVE_END = "exclusive_end"
E_END = "exclusive_end"
ION = "ion"
PROTEIN_ID = "protein_id"
PRODUCT_ION_TABLE = "product_ions"
PROTEIN_TABLE = "proteins"
SUBSEQ = "subsequence"

# Numeric constants
MAX_PEPTIDE_LEN = 50
MAX_KMER_LEN = 50
B_ION_AS_INT = 0
Y_ION_AS_INT = 1
ION_CHARGES_TO_CONSIDER = [1, 2]

# For simplicity
samples = [f"BMEM_AspN_Fxn{val}" for val in [4, 5, 6, 7, 8, 9]]
THOMAS_SAMPLES = [f"BMEM_AspN_Fxn{val}" for val in [4, 5, 6, 7, 8, 9]]


# ## Asserts
# assert file_exists(COMET_PARAMS) == True, "Comet executable not found!"
# assert file_exists(COMET_EXECUTABLE) == True, "Comet executable not found!"
