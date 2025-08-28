from pathlib import Path

from src.constants import COMET_DIR, DEFAULT_CRUX_PARAMS, MOUSE_PROTEOME
from src.crux import Crux, HSConfig


def test_create_native_run_comet_config(tmp_path, test_data_dir):

    name = "mouse_testing"
    mzml_to_scans = {
        # Path("data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML"): [0],
        # Path("data/spectra/mouse_samples/BMEM_AspN_Fxn5.mzML"): [0]
        test_data_dir / "spectra/10_mouse_spectra.mzML": [0],
        test_data_dir / "spectra/100_mouse_spectra.mzML": [0],
    }
    crux_comet_params = DEFAULT_CRUX_PARAMS
    fasta = MOUSE_PROTEOME
    # parent_out_dir = tmp_path
    parent_out_dir = Path("tmp")
    decoy_search = 2

    config = HSConfig(
        name=name,
        mzml_to_scans=mzml_to_scans,
        crux_comet_params=crux_comet_params,
        fasta=fasta,
        decoy_search=decoy_search,
        parent_out_dir=parent_out_dir,
        top_n_proteins=50,
    )
    config.create_native_run_comet_config()


def test_run_comet(tmp_path, test_data_dir):
    Crux(path=COMET_DIR / "crux-4.3.Darwin.x86_64/bin/crux").run_comet(
        mzml=test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML",
        fasta=test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta",
        crux_comet_params=test_data_dir / "crux.comet.params",
        decoy_search=2,
        out_dir=tmp_path,
    )
