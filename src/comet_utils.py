import os
import subprocess
from pathlib import Path
from typing import List, Union

import pandas as pd
from pydantic import BaseModel, field_validator

from src.erik_constants import SAMPLE
from src.erik_utils import make_dir


class CometExe(BaseModel):
    exe: Path
    params: Path


class CometResult(BaseModel):
    """Class for rows of Comet output"""

    sample: str
    scan: int
    ions_matched: int
    proteins: List[str]
    protein_count: int

    @field_validator("proteins", mode="before")
    def split_protein_column_by_comma(cls, protein: str) -> List[str]:
        return protein.split(",")

    # @model_validator(mode="after")
    # def matching_num_of_proteins(self):
    #     if len(self.proteins) != self.protein_count:
    #         print(
    #             f"Protein count = {self.protein_count} doesn't match the length of protein list = {len(self.proteins)}:\n{self.proteins} "
    #         )
    #     return self

    @classmethod
    def from_txt(cls, file_path: str) -> List["CometResult"]:
        file_path = Path(file_path)
        df = pd.read_csv(file_path, sep="\t", header=1)
        df["sample"] = file_path.stem
        return cls.from_dataframe(comet_output_df=df)

    @classmethod
    def from_dataframe(cls, comet_output_df: pd.DataFrame) -> List["CometResult"]:

        return [
            cls(
                sample=row["sample"],
                scan=row["scan"],
                ions_matched=row["ions_matched"],
                protein_count=row["protein_count"],
                proteins=row["protein"],
            )
            for _, row in comet_output_df.iterrows()
        ]


def read_comet_txt_to_df(txt_path: Path, sample: Union[None, str] = None):
    df = pd.read_csv(txt_path, sep="\t", header=1)
    # Add a "sample" column
    if sample is None:
        sample = txt_path.stem
    df[SAMPLE] = sample
    return df


# Run Comet
def run_comet_and_save_params_file(
    mzml_path: Path, comet: CometExe, parent_output_dir: Path
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
    make_dir(parent_output_dir)
    results_output_dir = parent_output_dir / mzml_path.stem
    make_dir(results_output_dir)

    # Copy the comet.params file the output directory
    subprocess.run(
        f"cp {comet.params} {results_output_dir}",
        shell=True,
        capture_output=True,
        text=True,
    )

    return run_comet(
        mzml_path=mzml_path,
        comet=comet,
        output_dir=results_output_dir,
        file_name_stem=mzml_path.stem,
    )


def run_comet(mzml_path: Path, comet: CometExe, output_dir: Path, file_name_stem: str):
    """
    Runs Comet on given MZML file and saves the results in the output_dir folder with
    file_name_stem as the result file names' stems. Typically there will be two Comet outputs:
        (1) output_dir/file_name_stem.txt
        (2) output_dir/file_name_stem.pep.xml
    """
    assert output_dir.exists(), f"Output dir {output_dir} doesn't exist!"

    # When running Comet, we need the CWD to contain the comet.params file
    os.chdir(comet.params.parent)
    output_path = output_dir / file_name_stem

    # Run Comet
    command = f"{comet.exe} {mzml_path} -N{output_path}"
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True,
    )

    # Check that expected output files exist -- TODO: this should be moved to a test
    expected_output_files = (
        Path(str(output_path) + ".txt"),
        Path(str(output_path) + ".pep.xml"),
    )
    for f in expected_output_files:
        assert f.exists(), f"{f} is expected to exist but does not"
    return result
