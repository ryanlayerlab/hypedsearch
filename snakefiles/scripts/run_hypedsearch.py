from pathlib import Path

from snakemake.script import snakemake

from src.hypedsearch import create_and_score_hybrids_for_spectrum
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta
from src.utils import load_json, setup_logger

# Constants in Snakemake's input and wildcards
mzml = Path(snakemake.input.mzml)
scan = int(snakemake.wildcards.scan)
spectrum = Spectrum.get_spectrum(mzml=mzml, scan=scan)
sample = snakemake.wildcards.sample

# Constants in Snakemake's config
crux_path = Path(snakemake.config.crux_path)
out_dir = Path(snakemake.config.out_dir)
fasta = Path(snakemake.config.fasta)
protein_name_to_seq_map = Fasta(path=fasta).protein_name_to_seq_map
crux_comet_params = Path(snakemake.config.crux_comet_params)
kmer_db = KmerDatabase(db_path=Path(snakemake.config.kmer_db))
kmer_to_protein_map = load_json(Path(snakemake.config.kmer_to_protein_map))
num_peaks = int(snakemake.config.num_peaks)
precursor_mz_ppm_tol = int(snakemake.config.precursor_mz_ppm_tol)
peak_to_ion_ppm_tol = int(snakemake.config.peak_to_ion_ppm_tol)
min_cluster_len = int(snakemake.config.min_cluster_len)
min_cluster_support = int(snakemake.config.min_cluster_support)
log_dir = snakemake.config.log_dir

# Setup logger
if len(log_dir) > 0:
    log_file = Path(log_dir) / f"{Mzml.get_mzml_name(mzml=mzml)}.{scan}.log"
    logger = setup_logger(str(log_file))
    logger.info(f"Logging to {log_file}")
else:
    log_file = ""
    logger = setup_logger()

# Run
create_and_score_hybrids_for_spectrum(
    spectrum=spectrum,
    kmer_db=kmer_db,
    protein_name_to_seq_map=protein_name_to_seq_map,
    kmer_to_proteins_map=kmer_to_protein_map,
    precursor_mz_ppm_tol=precursor_mz_ppm_tol,
    peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
    min_cluster_len=min_cluster_len,
    min_cluster_support=min_cluster_support,
    crux_comet_params=crux_comet_params,
    out_dir=out_dir,
    fasta=fasta,
    crux_path=crux_path,
    num_peaks=num_peaks,
)
if len(log_dir) > 0:
    logger.info("Deleteing log file...")
    log_file.unlink()
