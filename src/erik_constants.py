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
MOUSE_PROTEOME = GIT_REPO_DIR / "fastas/Uniprot_mouse.fasta"

## Strings
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

## Numeric constants
MAX_PEPTIDE_LEN = 50

## For simplicity
samples = [f"BMEM_AspN_Fxn{val}" for val in [4, 5, 6, 7, 8, 9]]


# ## Asserts
# assert file_exists(COMET_PARAMS) == True, "Comet executable not found!"
# assert file_exists(COMET_EXECUTABLE) == True, "Comet executable not found!"
