import logging
import os
import platform
import subprocess
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from time import time
from typing import Any, Dict, List, Literal, Optional, Union

import click
import pandas as pd
from pydantic import BaseModel, field_validator

from src.config import AppConfig
from src.constants import (
    COMET_DIR,
    COMET_PARAMS,
    COMET_RUN_1_DIR,
    COMET_RUN_2_DIR,
    DELTA_CN,
    EVAL,
    IONS_MATCHED,
    MOUSE_PROTEOME,
    NUM,
    PLAIN_PEPTIDE,
    PROTEIN,
    PROTEIN_COUNT,
    SAMPLE,
    SCAN,
    THOMAS_SAMPLES,
    XCORR,
)
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.peptides_and_ions import Peptide
from src.utils import (
    PathType,
    flatten_list_of_lists,
    log_params,
    make_directory,
    remove_gene_name,
    setup_logger,
)

MAC_OS = "macos"
LINUX_OS = "linux"

logger = logging.getLogger(__name__)


def get_comet_txts_in_dir(dir: Path) -> List[Path]:
    txt_files = list(dir.glob("*.txt"))
    pepxml_stems = {f.name[:-8] for f in dir.glob("*.pep.xml")}

    matching_txt_files = [
        txt_file for txt_file in txt_files if txt_file.stem in pepxml_stems
    ]

    return matching_txt_files


def read_comet_txt_to_df(txt_path: Path, sample: Optional[str] = None):
    """
    Reads Comet's output TXT file to a pandas dataframe.
    Adds a 'sample' column that's not in the TXT file but is in the XML file that Comet
    also spits out.
    """
    if sample is None:
        # Get sample name from the Comet XML file because
        xml_path = Path(str(txt_path).replace(".txt", ".pep.xml"))
        tree = ET.parse(xml_path)
        root = tree.getroot()
        namespace = {"pep": "http://regis-web.systemsbiology.net/pepXML"}
        msms_run_summary = root.find("pep:msms_run_summary", namespace)
        sample = msms_run_summary.get("base_name").split("/")[-1]

    # Read the TXT file to a dataframe and add the sample column
    df = pd.read_csv(txt_path, sep="\t", header=1)
    df[SAMPLE] = sample
    return df


class CometPSM(BaseModel):
    """Class for rows of Comet output"""

    sample: str
    num: int
    scan: int
    proposed_peptide: str
    ions_matched: int
    proteins: List[str]
    protein_count: int
    xcorr: float
    eval: float
    delta_cn: float
    # prev_aa: str
    # next_aa: str

    @field_validator("proteins", mode="before")
    def split_protein_column_by_comma(cls, protein: str) -> List[str]:
        """
        The 'protein' column in Comet's TXT output has string-valued values that look like:
            "<protein 1 name>,<protein 2 name>,<...>"
        I.e., comma-separated protein names.
        This function breaks the string of comma-separated values into a list:
            ["<protein 1 name>", "<protein 2 name>", ...]
        """
        return protein.split(",")

    @classmethod
    def from_txt(
        cls, file_path: str, as_df: bool = False
    ) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads Comet results .txt file to a list of dataclasses or a dataframe
        """
        file_path = Path(file_path)
        df = read_comet_txt_to_df(txt_path=file_path)
        if as_df:
            return df
        else:
            return cls.from_dataframe(comet_output_df=df)

    @classmethod
    def from_dataframe(cls, comet_output_df: pd.DataFrame) -> List["CometPSM"]:
        """
        Given a dataframe from the read_comet_txt_to_df method, convert the dataframe
        into a list of class instances--one for each row in the dataframe.
        """
        return [
            cls(
                sample=row[SAMPLE],
                scan=row[SCAN],
                num=row[NUM],
                ions_matched=row[IONS_MATCHED],
                protein_count=row[PROTEIN_COUNT],
                proteins=row[PROTEIN],
                proposed_peptide=row[PLAIN_PEPTIDE],
                xcorr=row[XCORR],
                eval=row[EVAL],
                delta_cn=row[DELTA_CN],
            )
            for _, row in comet_output_df.iterrows()
        ]

    def get_corresponding_spectrum(self) -> Spectrum:
        """
        Gets the spectrum corresponding to the Comet-returned PSM
        """
        return get_specific_spectrum_by_sample_and_scan_num(
            sample=self.sample, scan_num=self.scan
        )


def get_comet_protein_counts(
    comet_rows: Optional[List[CometPSM]] = None,
    shorten_names: bool = True,
) -> Counter:
    if comet_rows is None:
        comet_rows = load_comet_data(as_df=False)

    all_prots = flatten_list_of_lists([row.proteins for row in comet_rows])

    if shorten_names is True:
        all_prots = [remove_gene_name(protein_name=prot) for prot in all_prots]

    return Counter(all_prots)


def load_comet_data(
    samples: List[str] = THOMAS_SAMPLES, run: int = 1, as_df: bool = True
):
    """ """
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
        return CometPSM.from_dataframe(comet_output_df=comet_df)


@dataclass
class CometParams:
    path: Union[str, Path]
    file_lines: List[str] = field(init=False)

    def __post_init__(self):
        self.path = Path(self.path).absolute()
        assert self.path.exists(), "Path does not exist!"
        with open(self.path, "r") as f:
            self.file_lines = f.readlines()

    def update_database_name(self, fasta_path: Path):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("database_name ="):
                self.file_lines[line_idx] = f"database_name = {fasta_path}\n"
                break

    def update_scan_range(self, min_scan: int = 0, max_scan: int = 0):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("scan_range ="):
                self.file_lines[line_idx] = f"scan_range = {min_scan} {max_scan}\n"
                break

    def update_num_output_lines(self, num_output_lines: int):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("num_output_lines ="):
                self.file_lines[line_idx] = f"num_output_lines = {num_output_lines}\n"
                break

    def write(self, output_path: Path):
        with open(output_path, "w") as f:
            f.writelines(self.file_lines)


def get_os():
    os = platform.platform()
    if "macOS" in os:
        return MAC_OS
    elif "Linux" in os:
        return LINUX_OS
    else:
        raise RuntimeError(f"Unrecognized operating system: {os}")


def get_default_comet_executable_path():
    try:
        os = get_os()
        comet_path = COMET_DIR / f"comet.{os}.exe"
        assert comet_path.exists()
        return comet_path
    except:
        raise RuntimeError(
            "Please provide the Comet executable path. Could NOT find it in the default locations"
        )


@click.command(
    name="run-comet",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help=(
        "Path to .mzML file on which to run Comet. "
        "Or, if the path provided is to a directory, runs Comet on all the .mzML files found there"
    ),
)
@click.option(
    "--output_dir",
    "-o",
    type=PathType(),
    required=True,
    help="Path to directory where Comet outputs will be saved.",
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    default=MOUSE_PROTEOME,
    show_default=True,
    help="Path to FASTA file with which to run Comet",
)
@click.option(
    "--comet_exe",
    "-ce",
    type=PathType(),
    show_default=True,
    default=get_default_comet_executable_path(),
    help="Path to Comet executable",
)
@click.option(
    "--comet_params",
    "-cp",
    type=PathType(),
    show_default=True,
    default=COMET_PARAMS,
    help="Path to Comet parameters file, typically named comet.params",
)
@click.option(
    "--stem",
    "-n",
    type=str,
    help="The Comet output files will have this file stem. Defaults to the stem of the .mzML file",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    help=(
        "The scan number to run Comet on. "
        "Defaults to the scan_range in the provided Comet params file which is usually all scans in the .mzML file"
    ),
)
@click.option(
    "--num_psms",
    "-np",
    type=int,
    help=(
        "The number of PSMs per scan that Comet reports. "
        "Defaults to the num_output_lines in the provided Comet params file."
    ),
)
@click.option(
    "--check",
    "-c",
    is_flag=True,
    help="Check whether the Comet outputs already exist and, if they do, don't run Comet",
)
def cli(
    output_dir: Path,
    mzml: Path,
    fasta: Path = MOUSE_PROTEOME,
    stem: Optional[str] = None,
    comet_exe: Path = get_default_comet_executable_path(),
    comet_params: Path = COMET_PARAMS,
    check: bool = False,
    scan: Optional[int] = None,
    num_psms: Optional[int] = None,
):
    if mzml.is_dir():
        mzmls = list(mzml.glob("*.mzML"))
        logger.info(
            "You set the mzml argument to a directory. So running Comet on all .mzML files found there:\n"
            f"{mzmls}"
        )
        for mzml in mzmls:
            run_comet_on_one_mzml(
                mzml=mzml.absolute(),
                fasta=fasta.absolute(),
                comet_exe=comet_exe.absolute(),
                comet_params=comet_params.absolute(),
                output_dir=output_dir.absolute(),
                stem=stem,
                check=check,
                scan=scan,
                num_psms=num_psms,
            )

    else:
        run_comet_on_one_mzml(
            mzml=mzml.absolute(),
            fasta=fasta.absolute(),
            comet_exe=comet_exe.absolute(),
            comet_params=comet_params.absolute(),
            output_dir=output_dir.absolute(),
            stem=stem,
            check=check,
            scan=scan,
            num_psms=num_psms,
        )


def run_comet_on_one_mzml(
    output_dir: Path,
    mzml: Path,
    # keep_params: bool = True,
    fasta: Path = MOUSE_PROTEOME,
    stem: Optional[str] = None,
    comet_exe: Path = get_default_comet_executable_path(),
    comet_params: Path = COMET_PARAMS,
    check: bool = False,
    scan: Optional[int] = None,
    num_psms: Optional[int] = None,
) -> Path:
    start_time = time()
    # Constants & make sure paths are Path objects
    orig_dir = Path(os.getcwd()).absolute()
    comet_params_str = "comet.params"
    new_comet_params_path = (output_dir / comet_params_str).absolute()

    # Check whether expected output files already exist
    if stem is None:
        stem = f"{mzml.stem}"
    output_file_without_extension = (output_dir / f"{stem}").absolute()
    scan_range = ""
    if scan is not None:
        scan_range = f".{scan}-{scan}"
    comet_txt_output_path = Path(
        str(output_file_without_extension) + f"{scan_range}.txt"
    ).absolute()
    comet_pep_xml_output_path = Path(
        str(output_file_without_extension) + f"{scan_range}.pep.xml"
    ).absolute()
    if comet_txt_output_path.exists() and check:
        logger.info(
            f"Output file {comet_txt_output_path} already exists so Comet will NOT be run."
        )
        return comet_txt_output_path

    # Make sure parent_output_dir and sample_output_dir exist
    make_directory(output_dir)

    # Change to to sample_output_dir
    os.chdir(output_dir)

    # Copy template comet.params file to the sample_output_dir
    cmd = f"cp {comet_params} {new_comet_params_path}"
    cmd_result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    # Update comet.params with user-provided values
    params = CometParams(path=new_comet_params_path)
    params.update_database_name(fasta_path=fasta)
    if scan is not None:
        params.update_scan_range(min_scan=scan, max_scan=scan)
    if num_psms is not None:
        params.update_num_output_lines(num_output_lines=num_psms)
    params.write(output_path=new_comet_params_path)

    # Run Comet
    cmd = f"{comet_exe} {mzml} -N{output_file_without_extension}"
    logger.info(f"Running Comet with command:\n{cmd}\nfrom directory:\n{os.getcwd()}")
    cmd_result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    # Change back to original directory
    os.chdir(orig_dir)

    # Check for errors
    if cmd_result.returncode != 0:
        raise RuntimeError(
            f"Comet failed with stdout:\n{cmd_result.stdout}\n"
            f"And stderr:\n{cmd_result.stderr}"
        )

    # Check that expected output files exist -- TODO: this should be moved to a test
    expected_output_files = [comet_txt_output_path, comet_pep_xml_output_path]
    for f in expected_output_files:
        assert f.exists(), f"{f} is expected to exist but does not"
    logger.info(f"Running Comet took {round(time() - start_time, 2)} seconds")
    return comet_txt_output_path


if __name__ == "__main__":
    setup_logger()
    cli()
