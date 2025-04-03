import shutil
from pathlib import Path

import pytest

from hypedsearch import HypedsearchConfig, hypedsearch
from src.constants import COMET_PARAMS
from src.utils import file_hash, make_directory


def test_hypedsearch(tmp_path):
    parent_output_dir = Path("tmp/test").absolute()
    tmp_spectra_dir = parent_output_dir / "spectra"
    tmp_output_dir = parent_output_dir / "output"
    tmp_db_path = tmp_output_dir / "db.db"
    make_directory(tmp_spectra_dir)
    make_directory(tmp_output_dir)
    mzml_dir = Path("data/spectra").absolute()

    mzml_paths = list(mzml_dir.glob("*.mzML"))
    # mzml_paths = mzml_paths[:2]
    mzml_path = mzml_paths[0]
    for mzml in mzml_paths:
        shutil.copy(
            src=mzml,
            dst=tmp_spectra_dir,
        )
    config = HypedsearchConfig(
        mzml_dir=tmp_spectra_dir,
        output_dir=tmp_output_dir,
        top_n_proteins=50,
        db_path=tmp_db_path,
        mzml_path=mzml_path,
        scan_num=1,
    )

    hypedsearch(config)
