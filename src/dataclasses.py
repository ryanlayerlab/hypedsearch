import logging
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from itertools import groupby
from pathlib import Path
from typing import (
    Callable,
    Counter,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    Self,
    Set,
    Tuple,
    Union,
)

import pandas as pd
import seaborn as sns
from pydantic import BaseModel
from scipy.stats import ecdf, percentileofscore

from src.comet_utils import CometPSM, HybridPeptide
from src.constants import HYBRID, NATIVE
from src.kmer_database import KmerToProteinsMap
from src.peptides_and_ions import Peptide
from src.utils import flatten_list_of_lists

logger = logging.getLogger(__name__)


@dataclass
class CometPSMs:
    psms: List[CometPSM]

    def __len__(self):
        return len(self.psms)

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "Self":
        return cls(psms=CometPSM.from_txt(txt=path))

    def __post_init__(self):
        if isinstance(self.psms, (str, Path)):
            self.psms = CometPSM.from_txt(txt=self.psms)

    def get_psms_for_sample_and_scan(self, sample: str, scan: int):
        return list(
            filter(lambda psm: (psm.sample == sample) and (psm.scan == scan), self.psms)
        )

    def get_high_confidence_psms(self, q_value_threshold: float) -> "Self":
        return CometPSMs(
            psms=[psm for psm in self.psms if psm.q_value <= q_value_threshold]
        )

    def filter_out_non_hybrids(self) -> "Self":
        return CometPSMs(psms=[psm for psm in self.psms if psm.is_hybrid])

    @property
    def scans(self) -> Set[Tuple[str, int]]:
        return {(psm.sample, psm.scan) for psm in self.psms}

    def get_psms_for_scans(self, scans: Iterable[Tuple[str, int]]) -> "Self":
        """
        Args:
            scans (Iterable[Tuple[str, int]]): Iterable of (sample, scan) tuples
        """
        return CometPSMs(
            psms=list(filter(lambda psm: (psm.sample, psm.scan) in scans, self.psms))
        )

    @cached_property
    def protein_counts(self) -> Counter:
        all_comet_proteins = flatten_list_of_lists([psm.proteins for psm in self.psms])
        comet_protein_counts = Counter(all_comet_proteins)
        return comet_protein_counts

    @cached_property
    def kmer_counts(self) -> Dict[int, Counter]:
        kmer_counter = defaultdict(list)
        for psm in self.psms:
            for kmer in Peptide(seq=psm.seq).kmers(min_k=1, max_k=len(psm.seq)):
                kmer_counter[len(kmer.seq)].append(kmer.seq)
        return {k: Counter(seqs) for k, seqs in kmer_counter.items()}


@dataclass
class PSMScoring:
    name: str
    values_by_charge: Dict[int, List[float]] = field(
        default_factory=lambda: defaultdict(list)
    )

    def add_value(self, value: float, charge: int):
        self.values_by_charge[charge].append(value)

    def set_ecdfs(self):
        ecdfs_by_charge = dict()
        for charge, values in self.values_by_charge.items():
            ecdfs_by_charge[charge] = ecdf(values)
        self.ecdfs_by_charge = ecdfs_by_charge
