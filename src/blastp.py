import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

from src.peptides_and_ions import Fasta, Peptide


def run_blastp(
    query_peptides: Union[List[Peptide], List[str]],
    fasta: Union[str, Path],
    out_path: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    # Columns to include in BLAST output
    # See here for BLAST column help: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    blast_columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qcovhsp",
    ]
    # blast_columns = " ".join(blast_columns)

    # Set query peptides depending on input type
    if isinstance(query_peptides[0], str):
        query_peptides = [Peptide(seq=seq, name=seq) for seq in query_peptides]

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        # Create a FASTA consisting of query peptides
        peptide_fasta = tmp_dir / "peptides.faa"
        Fasta.write_fasta(peptides=query_peptides, out_path=peptide_fasta)
        # Run blastp
        if out_path is None:
            out_path = tmp_dir / "peptides.csv"
        cmd_parts = [
            "blastp",
            f"-query {peptide_fasta}",
            f"-db {fasta}",
            f"-out {out_path}",
            f'-outfmt "10 {" ".join(blast_columns)}"',
        ]
        cmd = " ".join(cmd_parts)
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            shell=True,
        )
        assert result.returncode == 0, f"blastp failed with stderr:\n{result.stderr}"

        # Load results
        df = pd.read_csv(out_path, names=blast_columns)
    return df
