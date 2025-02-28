from enum import Enum
from pathlib import Path

# from src.erik_utils import file_exists

# Paths
GIT_REPO_DIR = Path(__file__).parents[1]
ROOT_DIR = Path(__file__).parents[2]
TEST_DIR = GIT_REPO_DIR / "tests"
LOGS_DIR = GIT_REPO_DIR / "logs"
PLOTS_DIR = GIT_REPO_DIR / "plots"
DATA_DIR = GIT_REPO_DIR / "data"
RESULTS_DIR = DATA_DIR / "results"
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
NEUTRAL_MASS = "neutral_mass"
KMERS = "kmers"
SEQ = "sequence"
COUNT = "count"
BOTH = "both"
KMER_AND_MASS = "kmer+mass"
KMER = "kmer"
SAMPLE = "sample"
IONS_MATCHED = "ions_matched"
SPECTRUM_ID = "spectrum_id"
PROTEIN_COUNT = "protein_count"
PROTEIN = "protein"
XCORR = "xcorr"
EVAL = "e-value"
PROPOSED_PROTEIN = "proposed_protein"
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
B_ION_TYPE = "b"
Y_ION_TYPE = "y"

MEMORY = ":memory:"


class IonTypes(Enum):
    B_ION_TYPE = "b"
    Y_ION_TYPE = "y"


ION_TYPE_TO_INT = {
    "b": 0,
    "y": 1,
}

ION_INT_TO_TYPE = {v: k for k, v in ION_TYPE_TO_INT.items()}

ALL_IONS = "all"

# Numeric constants
MAX_PEPTIDE_LEN = 50
MAX_KMER_LEN = 50
DEFAULT_MAX_K = 25
DEFAULT_MIN_K = 1
B_ION_AS_INT = 0
Y_ION_AS_INT = 1
ION_CHARGES_TO_CONSIDER = [1, 2]

# For simplicity
samples = [f"BMEM_AspN_Fxn{val}" for val in [4, 5, 6, 7, 8, 9]]
THOMAS_SAMPLES = [f"BMEM_AspN_Fxn{val}" for val in [4, 5, 6, 7, 8, 9]]


# ## Asserts
# assert file_exists(COMET_PARAMS) == True, "Comet executable not found!"
# assert file_exists(COMET_EXECUTABLE) == True, "Comet executable not found!"

# Amino acid masses
# Source for dictionary below is the code here https://www.rapidnovor.com/mxie-customized/massCal.js
# which is used by this website to compute masses: https://www.rapidnovor.com/mass-calculator/
AMINO_ACID_MASSES = {
    "A": 71.0371138,
    "R": 156.1011110,
    "N": 114.0429274,
    "D": 115.0269430,
    "C": 103.0091845,
    "E": 129.0425931,
    "Q": 128.0585775,
    "G": 57.0214637,
    "H": 137.0589119,
    "I": 113.0840640,
    "L": 113.0840640,
    "K": 128.0949630,
    "M": 131.0404846,
    "F": 147.0684139,
    "P": 97.0527639,
    "S": 87.0320284,
    "T": 101.0476785,
    "W": 186.0793130,
    "Y": 163.0633285,
    "U": 168.9641990,
    "V": 99.0684139,
}

# Other masses
HYDROGEN_MASS = 1.007825035
WATER_MASS = 18.010564686  # from https://www.rapidnovor.com/mxie-customized/massCal.js
OXYGEN_MASS = 15.99491463
PROTON_MASS = 1.00727646688  # from Scott
# This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
WATER_MASS = (2 * HYDROGEN_MASS) + OXYGEN_MASS
