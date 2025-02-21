import re
from pathlib import Path
from typing import List, Optional, Union

from pydantic import BaseModel


def to_path(path: Union[str, Path]):
    return Path(path)


class Position(BaseModel):
    inclusive_start: Optional[int] = None
    exclusive_end: Optional[int] = None


class Kmer(BaseModel):
    seq: str
    position: Position


def generate_aa_kmers(aa_seq: str, max_k: int, min_k: Optional[int] = 1) -> List[Kmer]:
    disallowed_amino_acid_symbols = ["B", "X", "U", "Z", "O", "J"]
    kmers = []
    for kmer_len in range(min_k, max_k + 1):  # Iterate over all k from 1 to max_k
        for inclusive_start in range(
            len(aa_seq) - kmer_len + 1
        ):  # Slide over the string
            exclusive_end = inclusive_start + kmer_len
            kmer = aa_seq[inclusive_start:exclusive_end]
            if any(char in disallowed_amino_acid_symbols for char in kmer):
                continue
            kmers.append(
                Kmer(
                    seq=kmer,
                    position=Position(
                        inclusive_start=inclusive_start,
                        exclusive_end=exclusive_end,
                    ),
                )
            )

    return kmers


def get_positions_of_subseq_in_seq(subseq: str, seq: str) -> List[Position]:
    matches = re.finditer(subseq, seq)
    return [
        Position(inclusive_start=match.start(), exclusive_end=match.end())
        for match in matches
    ]
