from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from src.comet_utils import CometPSM, get_comet_protein_counts, load_comet_data
from src.constants import (
    COMET_COUNTS,
    COMET_RUN_1_DIR,
    FULL_NAME,
    MOUSE_PROTEOME,
    PROTEIN,
    QUANTILE,
    RESULTS_DIR,
    SHORT_NAME,
    SPECTRA_DIR,
    VALIDATED,
)
from src.peptides_and_ions import Peptide
from src.utils import flatten_list_of_lists, pickle_and_compress, remove_gene_name


@dataclass
class SampleFiles:
    validated_proteins: Optional[Path] = None
    mzml: Optional[Path] = None
    comet_run_1: Optional[Path] = None
    comet_run_2: Optional[Path] = None


def get_sample_files(sample_num: int):
    recognized_samples = list(range(4, 10))
    assert sample_num in recognized_samples

    comet_run_1 = (
        COMET_RUN_1_DIR / f"BMEM_AspN_Fxn{sample_num}/BMEM_AspN_Fxn{sample_num}.txt"
    )
    validated = (
        SPECTRA_DIR
        / f"rehybridpeptidelist/Specmill_Validated_Proteins_Fxn0{sample_num}.txt"
    )
    mzml = SPECTRA_DIR / f"BMEM_AspN_Fxn{sample_num}.mzML"
    return SampleFiles(
        validated_proteins=validated,
        comet_run_1=comet_run_1,
        mzml=mzml,
    )


def get_protein_name_from_validated_proteins_row(row: pd.Series) -> str:
    """
    The input is a row of Thomas's validated protein 'Specmill_Validated_Proteins_<sample>.txt' file.
    Returns the name of the protein corresponding to the row in this form:
        <row.database>|<row.accession_number>.

    Only the SwissProt=sp database is supported because that's all that's been in Thomas's
    files so far
    """
    recognized_db_names = {"SwissProt": "sp"}
    database = row["database"]
    accession = row["accession_number"]

    return f"{recognized_db_names[database]}|{accession}"


def read_validated_proteins_txt(file_path: str, sample_num: int) -> pd.DataFrame:
    """
    Read one of Thomas's validated protein 'Specmill_Validated_Proteins_<sample>.txt' files
    to a dataframe
    """
    validated_proteins = pd.read_csv(file_path, sep="\t")

    # Remove "enzyme" column which seems empty so it's being read incorrectly
    if sample_num != 8:
        enzyme_idx = list(validated_proteins.columns).index("enzyme")
        columns = list(validated_proteins.columns)
        _ = columns.pop(enzyme_idx)
        validated_proteins = pd.DataFrame(
            data=validated_proteins.iloc[:, :-1].to_numpy(), columns=columns
        )

    # Add protein name column to dataframe: <database>|<accession number>
    validated_proteins[SHORT_NAME] = validated_proteins.apply(
        lambda row: get_protein_name_from_validated_proteins_row(row), axis=1
    )

    return validated_proteins


def get_all_validated_proteins() -> pd.DataFrame:
    valid_proteins_dfs = []
    for sample_num in range(4, 10):
        files = get_sample_files(sample_num=sample_num)
        validated_proteins = read_validated_proteins_txt(
            file_path=files.validated_proteins, sample_num=sample_num
        )
        valid_proteins_dfs.append(validated_proteins)
    valid_prots_df = pd.concat(valid_proteins_dfs, ignore_index=True)
    return valid_prots_df


def comet_counts_and_valid_prots(fasta_path: str = MOUSE_PROTEOME) -> pd.DataFrame:
    fasta_prots = set([prot.name for prot in Peptide.from_fasta(fasta_path=fasta_path)])
    shortened_fasta_prots = set([remove_gene_name(prot) for prot in fasta_prots])
    valid_prots = set(get_all_validated_proteins()[SHORT_NAME])
    print(
        (
            f"Here are the proteins that are in the validated list but NOT in the FASTA: {valid_prots.difference(shortened_fasta_prots)}"
        )
    )
    comet_prot_counts = get_comet_protein_counts(shorten_names=False)
    df = defaultdict(list)
    for prot in fasta_prots:
        short_name = remove_gene_name(protein_name=prot)
        df[FULL_NAME].append(prot)
        df[SHORT_NAME].append(short_name)
        df[COMET_COUNTS].append(comet_prot_counts[prot])
        df[VALIDATED].append(int(short_name in valid_prots))
    return pd.DataFrame(df)
