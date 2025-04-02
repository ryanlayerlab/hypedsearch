import logging
import os
import subprocess
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from time import time
from typing import List, Optional, Union

import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import (
    COMET_RUN_1_DIR,
    COMET_RUN_2_DIR,
    DELTA_CN,
    EVAL,
    IONS_MATCHED,
    MOUSE_PROTEOME,
    PLAIN_PEPTIDE,
    PROTEIN,
    PROTEIN_COUNT,
    SAMPLE,
    SCAN,
    THOMAS_SAMPLES,
    XCORR,
)
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.peptide_spectrum_comparison import PeptideSpectrumComparison
from src.peptides_and_ions import Peptide
from src.utils import (
    flatten_list_of_lists,
    log_params,
    make_directory,
    remove_gene_name,
)

logger = logging.getLogger(__name__)


def read_comet_txt_to_df(txt_path: Path):
    """
    Reads Comet's output TXT file to a pandas dataframe.
    Adds a 'sample' column that's not in the TXT file but is in the XML file that Comet
    also spits out.
    """
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

    def compare(self, spectrum: Optional[Spectrum] = None):
        if spectrum is None:
            logger.info(
                "Spectrum not provided. "
                "So getting the spectrum which might be a problem since no peak filtering will happen."
            )
            spectrum = self.get_corresponding_spectrum()
        psm = PeptideSpectrumComparison(
            spectrum=spectrum, peptide=Peptide(seq=self.proposed_peptide)
        )
        psm.compare()
        return psm


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


def update_database_path_in_comet_params_file(
    fasta_path: str, comet_params_path: str, output_path: str
) -> None:
    """
    Given the path to a comet.params file, update the "database_name=..." line
    to point to the given FASTA path instead of what it was set to. Save the resulting,
    updated comet.params file to the given output path.
    """
    with open(comet_params_path, "r") as infile:
        lines = infile.readlines()

    # Find and replace the line that starts with "database_name ="
    for line_idx, line in enumerate(lines):
        if line.strip().startswith("database_name ="):
            lines[line_idx] = f"database_name = {fasta_path}\n"
            break  # Stop after finding the first occurrence

    # Write the modified content to the output path
    with open(output_path, "w") as outfile:
        outfile.writelines(lines)


@log_params
def run_comet(
    template_comet_params_path: Union[str, Path],
    fasta_path: Union[str, Path],
    mzml_path: Union[str, Path],
    output_dir: Union[str, Path],
    comet_exe_path: Union[str, Path],
    keep_params: bool = True,
    overwrite: bool = True,
    output_file_stem: Optional[str] = None,
) -> Path:
    start_time = time()
    # Constants & make sure paths are Path objects
    orig_dir = Path(os.getcwd()).absolute()
    try:
        comet_params_str = "comet.params"
        sample = mzml_path.stem
        mzml_path = Path(mzml_path).absolute()
        output_dir = Path(output_dir).absolute()
        comet_exe_path = Path(comet_exe_path).absolute()
        comet_params_path = (output_dir / comet_params_str).absolute()

        # Check whether expected output files already exist
        if output_file_stem is None:
            output_file_stem = f"{sample}"
        output_file_without_extension = (output_dir / f"{output_file_stem}").absolute()
        comet_txt_output_path = Path(
            str(output_file_without_extension) + ".txt"
        ).absolute()
        comet_pep_xml_output_path = Path(
            str(output_file_without_extension) + ".pep.xml"
        ).absolute()
        if comet_txt_output_path.exists() and not overwrite:
            logger.info(
                f"Output file {comet_txt_output_path} already exists and overwrite is set to {overwrite}. "
                f"Returning early."
            )
            return comet_txt_output_path

        # Make sure parent_output_dir and sample_output_dir exist
        make_directory(output_dir)

        # Change to to sample_output_dir
        os.chdir(output_dir)

        # Copy template comet.params file to the sample_output_dir
        cmd = f"cp {template_comet_params_path} {comet_params_path}"
        cmd_result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
        )

        # Update "database_name..." line in comet.params to point to given FASTA file
        update_database_path_in_comet_params_file(
            fasta_path=fasta_path,
            comet_params_path=comet_params_path.absolute(),
            output_path=comet_params_path,
        )

        # Run Comet
        cmd = f"{comet_exe_path} {mzml_path} -N{output_file_without_extension}"
        logger.info(
            f"Running Comet with command:\n{cmd}\nfrom directory:\n{os.getcwd()}"
        )
        cmd_result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
        )

        # Check that expected output files exist -- TODO: this should be moved to a test
        expected_output_files = [comet_txt_output_path, comet_pep_xml_output_path]
        for f in expected_output_files:
            assert f.exists(), f"{f} is expected to exist but does not"

        logger.info(
            f"Looks like Comet ran successfully. The expected files ({expected_output_files}) exist"
        )

        # Delete comet.params if desired
        if not keep_params:
            comet_params_path.unlink()

        # Change back to original directory
        os.chdir(orig_dir)

        logger.info(f"Running Comet took {round(time() - start_time, 2)} seconds")
        return comet_txt_output_path

    except:
        os.chdir(orig_dir)
        raise RuntimeError(
            f"Running Comet failed. The command:\n{cmd}\nfailed with error:\n{cmd_result.stderr}"
        )
