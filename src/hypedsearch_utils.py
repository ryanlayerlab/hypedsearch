from pathlib import Path
from typing import List

import pandas as pd
from pydantic import BaseModel, field_validator


class HSResult(BaseModel):
    """Class for rows of HypedSearch output"""

    sample: str
    spectrum_id: int
    hybrid: bool  # True = hybrid, False = native
    seq: str
    precursor_mass: float
    left_kmers: List[str]
    right_kmers: List[str]

    @field_validator("left_kmers", "right_kmers", mode="before")
    @classmethod
    def split_by_backslash(cls, field_val: str) -> List[str]:
        """
        HypedSearch can return multiple proteins in '/' separated list in the
        right_kmer and left_kmer columns
        """
        return field_val.split("/")

    @classmethod
    def from_txt(cls, file_path: str) -> List["HSResult"]:
        file_path = Path(file_path)
        df = pd.read_csv(file_path, sep="\t")
        df["sample"] = file_path.stem[3:]  # ignore the "HS_" prefix
        return cls.from_dataframe(hypedsearch_output_df=df)

    @classmethod
    def from_dataframe(cls, hypedsearch_output_df: pd.DataFrame) -> List["HSResult"]:

        return [
            cls(
                sample=row["sample"],
                spectrum_id=row["spectrum_id"],
                hybrid=(True if (row["hybrid"] == "Hybrid") else False),
                seq=row["sequence"],
                precursor_mass=row["precursor mass"],
                left_kmers=row["left kmer"],
                right_kmers=row["right kmer"],
            )
            for _, row in hypedsearch_output_df.iterrows()
        ]
