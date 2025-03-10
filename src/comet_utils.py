import logging
import os
import subprocess
from collections import Counter
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import (
    COMET_RUN_1_DIR,
    COMET_RUN_2_DIR,
    IONS_MATCHED,
    MOUSE_PROTEOME,
    PLAIN_PEPTIDE,
    PROTEIN,
    PROTEIN_COUNT,
    SAMPLE,
    SCAN,
    THOMAS_SAMPLES,
)
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.utils import (
    flatten_list_of_lists,
    log_params,
    make_directory,
    remove_gene_name,
)

logger = logging.getLogger(__name__)


class CometExe(BaseModel):
    exe: Path
    params: Path


class CometRow(BaseModel):
    """Class for rows of Comet output"""

    sample: str
    scan: int
    ions_matched: int
    proteins: List[str]
    protein_count: int
    proposed_peptide: str

    @field_validator("proteins", mode="before")
    def split_protein_column_by_comma(cls, protein: str) -> List[str]:
        return protein.split(",")

    @classmethod
    def from_txt(cls, file_path: str) -> List["CometRow"]:
        file_path = Path(file_path)
        df = pd.read_csv(file_path, sep="\t", header=1)
        df["sample"] = file_path.stem
        return cls.from_dataframe(comet_output_df=df)

    @classmethod
    def from_dataframe(cls, comet_output_df: pd.DataFrame) -> List["CometRow"]:

        return [
            cls(
                sample=row[SAMPLE],
                scan=row[SCAN],
                ions_matched=row[IONS_MATCHED],
                protein_count=row[PROTEIN_COUNT],
                proteins=row[PROTEIN],
                proposed_peptide=row[PLAIN_PEPTIDE],
            )
            for _, row in comet_output_df.iterrows()
        ]

    def get_corresponding_spectrum(self) -> Spectrum:
        return get_specific_spectrum_by_sample_and_scan_num(
            sample=self.sample, scan_num=self.scan
        )


def get_comet_protein_counts(
    comet_rows: Optional[List[CometRow]] = None, shorten_names: bool = True
) -> Counter:
    if comet_rows is None:
        comet_rows = load_comet_data(as_df=False)

    all_prots = flatten_list_of_lists([row.proteins for row in comet_rows])

    if shorten_names is True:
        all_prots = [remove_gene_name(protein_name=prot) for prot in all_prots]

    return Counter(all_prots)


def read_comet_txt_to_df(txt_path: Path, sample: Union[None, str] = None):
    df = pd.read_csv(txt_path, sep="\t", header=1)
    # Add a "sample" column
    if sample is None:
        sample = txt_path.stem
    df[SAMPLE] = sample
    return df


def run_comet_and_save_params_file(
    mzml_path: Path,
    comet: CometExe,
    parent_output_dir: Path,
    fasta_path: Path = MOUSE_PROTEOME,
):
    """
    Given an MZML path/name.mzML and an output dir (parent_output_dir), this function:
        1. creates a folder parent_output_dir/name if it doesn't already exist
        2. copies comet.params to the parent_output_dir/name directory; the reason to do
            this is to keep Comet results with the comet.params file that generated them
        3. runs Comet and saves the result files in parent_output_dir/name/name.<ext>.
            Typically there will be two Comet outputs:
            (1) parent_output_dir/name/name.txt
            (2) parent_output_dir/name/name.pep.xml
    """
    # Make sure the parent directory exists and make sure the
    make_directory(parent_output_dir)
    results_output_dir = parent_output_dir / mzml_path.stem
    make_directory(results_output_dir)

    # Copy the comet.params file the output directory
    subprocess.run(
        f"cp {comet.params} {fasta_path} {results_output_dir}",
        shell=True,
        capture_output=True,
        text=True,
    )

    exit_code = run_comet(
        mzml_path=mzml_path.absolute(),
        comet=comet,
        output_dir=results_output_dir,
        file_name_stem=mzml_path.stem,
    )

    if exit_code == 0:
        copied_fasta = results_output_dir / fasta_path.name
        logger.info(f"Deleting {copied_fasta}")
        copied_fasta.unlink()


def run_comet(mzml_path: Path, comet: CometExe, output_dir: Path, file_name_stem: str):
    """
    Runs Comet on given MZML file and saves the results in the output_dir folder with
    file_name_stem as the result file names' stems. Typically there will be two Comet outputs:
        (1) output_dir/file_name_stem.txt
        (2) output_dir/file_name_stem.pep.xml
    """
    assert output_dir.exists(), f"Output dir {output_dir} doesn't exist!"
    curr_dir = os.getcwd()

    # When running Comet, we need the CWD to contain the comet.params file
    os.chdir(output_dir.absolute())
    output_path = output_dir / file_name_stem

    # Run Comet
    command = f"{comet.exe} {mzml_path.absolute()} -N{output_path.absolute()}"
    logger.info(
        f"Running Comet with command:\n{command}\nfrom directory:\n{os.getcwd()}"
    )
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        # Comet failed so clean up and exit
        msg = (
            "Comet failed!\n"
            f"Exit code = {result.returncode}\n"
            f"Error output = \n{result.stderr}"
        )
        logger.info(msg)
        os.chdir(curr_dir)
        return result.returncode

    # Check that expected output files exist -- TODO: this should be moved to a test
    expected_output_files = (
        Path(str(output_path) + ".txt"),
        Path(str(output_path) + ".pep.xml"),
    )
    for f in expected_output_files:
        assert f.exists(), f"{f} is expected to exist but does not"

    logger.info(
        f"Looks like Comet ran successfully. The expected files ({expected_output_files}) exist"
    )
    os.chdir(curr_dir)
    return result.returncode


def load_comet_data(
    samples: List[str] = THOMAS_SAMPLES, run: int = 1, as_df: bool = True
):
    if run == 1:
        comet_results_dir = COMET_RUN_1_DIR
    elif run == 2:
        comet_results_dir = COMET_RUN_2_DIR
    else:
        raise ValueError(f"'run' must be 1 or 2. You set run={run}")
    comet_dfs = []
    for sample in samples:
        logger.info(f"Reading data for {sample}")
        comet_output = comet_results_dir / f"{sample}/{sample}.txt"
        comet_dfs.append(read_comet_txt_to_df(txt_path=comet_output))
    comet_df = pd.concat(comet_dfs, ignore_index=True)

    if as_df is True:
        return comet_df
    else:
        return CometRow.from_dataframe(comet_output_df=comet_df)
