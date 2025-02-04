from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class Kmer:
    seq: str
    inclusive_start: Optional[int] = None
    exclusive_end: Optional[int] = None


@dataclass
class KmerWithMass:
    kmer: Kmer
    mass: float


@dataclass
class Protein:
    seq: str
    id: Optional[int] = None
    desc: Optional[str] = None


@dataclass
class Peak:
    mz: float
    abundance: float


@dataclass
class Spectrum:
    # mass_over_charges: Tuple[float, ...]
    # abundances: Tuple[float, ...]
    peaks: Tuple[Peak, ...]
    precursor_mz: float
    precursor_charge: int
    precursor_abundance: float
    spectrum_id: str
    retention_time: float
    mzml: Optional[str] = None
    scan_num: Optional[int] = None


@dataclass(frozen=True)
class KmerTableRow:
    mass: float
    start: int
    end: int
    ion: int
    charge: int
    protein_id: int


@dataclass
class KmerIons:
    b_ion: KmerTableRow
    y_ion: KmerTableRow
