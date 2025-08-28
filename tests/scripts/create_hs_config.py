from src.hypedsearch import HybridRunConfig

mzml_to_scans = {
    "tests/data/spectra/10_mouse_spectra.mzML": [1, 2, 3],
    "tests/data/spectra/BMEM_AspN_Fxn4_scans1-20.mzML": [7, 10],
}
config = HybridRunConfig(
    mzml_to_scans=mzml_to_scans,
    out_dir="tests/data/create_db_tutorial/scan_results",
    fasta="fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta",
    kmer_db="tests/data/mouse_top_10_proteins.db",
    kmer_to_protein_lsmap="tests/data/create_db_tutorial/kmer_to_prot.json",
    crux_path="/usr/local/bin/crux",
)
config_path = "tests/data/create_db_tutorial/hs.snakemake.config"
print(f"Saving config to {config_path}")
config.save("tests/data/create_db_tutorial/hs.snakemake.config")
