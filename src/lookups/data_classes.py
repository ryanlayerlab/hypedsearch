from dataclasses import dataclass
from typing import Literal, Optional, Tuple


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
    sequence: str
    protein_id: Optional[int] = None
    desc: Optional[str] = None


@dataclass
class Peak:
    mz: float
    abundance: float
    id: Optional[int] = None


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
class ProductIonTableRow:
    mass: float
    start: int
    end: int
    ion: int
    charge: int
    protein_id: int


@dataclass(frozen=True)
class IonWithProteinInfo(ProductIonTableRow):
    subsequence: str


@dataclass
class KmerIons:
    b_ion: ProductIonTableRow
    y_ion: ProductIonTableRow


@dataclass(frozen=True)
class Ion:
    seq: str
    charge: int
    mz: float
    ion_type: Literal["b", "y"]


@dataclass(frozen=True)
class PeakIonMatch:
    peak: Peak
    ion: Ion
