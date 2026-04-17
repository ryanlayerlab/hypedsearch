"""
Classes and methods for working with peptide-spectrum matches (PSMs).
"""

import logging
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Union

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from pydantic import BaseModel
from scipy.interpolate import PchipInterpolator

from src.constants import (
    B_ION_TYPE,
    CALC_NEUTRAL_MASS,
    CHARGE,
    COMET,
    COMET_PROTEIN_SEPARATOR,
    CRUX,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_THRESHOLD,
    DELTA_CN,
    EVAL,
    EXP_NEUTRAL_MASS,
    IONS_MATCHED,
    IONS_TOTAL,
    NUM,
    PEPTIDE_NEUTRAL_MASS,
    PLAIN_PEPTIDE,
    PROTEIN,
    Q_VAL,
    RETENTION_TIME,
    RETENTION_TIME_STR,
    SAMPLE,
    SCAN,
    SPECTRUM_NEUTRAL_MASS,
    SPECTRUM_PRECURSOR_MZ,
    XCORR,
    Y_ION_TYPE,
)
from src.hybrids_via_clusters import HybridPeptide
from src.mass_spectra import Mzml, Peak, Spectrum, organize_by_spectrum_uid, plot_peaks
from src.peptides_and_ions import Fasta, Peptide, compute_peptide_precursor_mz
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    save_fig,
    set_title_axes_labels,
)
from src.utils import (
    flatten_list_of_lists,
    list_to_df,
    load_json,
    mass_difference_in_ppm,
    to_json,
)

logger = logging.getLogger(__name__)


@dataclass
class CometTxt:
    """
    Class for handling Comet output .txt files where from Comet directly or from `crux comet`.
    This class can handle both formats.
    """

    path: Union[str, Path]
    sample: str = field(init=False)
    file_type: Literal["comet", "crux"] = field(init=False)

    def __post_init__(self):
        # Set file_type
        with open(self.path, "r") as f:
            first_line = f.readline().strip()
            if ("comet" in first_line) or ("Comet" in first_line):
                self.file_type = COMET
            else:
                self.file_type = CRUX

        # Set sample
        self.sample = self.path.stem.split(".")[0]


class PeakIonMatch(BaseModel):
    """
    Class to represent a match between a theoretical product ion and a peak in a spectrum.
    """

    ion_mz: float
    ion_charge: int
    ion_seq: str
    ion_type: str
    peak_mz: float
    peak_intensity: float
    sample: Optional[str]
    scan: Optional[int]

    def mz_diff(self, type: Literal["rel", "rel_ppm"] = "rel_ppm") -> float:
        """
        Returns the m/z difference between the theoretical ion and the spectrum peak.
        Let x_i = theoretical ion mass, x_p = peak mass,
        then returns
            - (x_i - x_t) / x_i when type='rel'
            - ((x_i - x_t) / x_i) * (10**6) when type='rel_ppm'

        Args:
            type: Specifies whether to return the relative m/z different in PPM or not. Defaults to "rel_ppm".
        """
        if type == "rel":
            return (self.ion_mz - self.peak_mz) / self.ion_mz
        elif type == "rel_ppm":
            return ((self.ion_mz - self.peak_mz) / self.ion_mz) * (10**6)


def get_peaks_near_mz(
    query_mz: float, peaks: List[Peak], ppm_tolerance: float
) -> List[Peak]:
    """
    Given a list of mass spectrum peaks and a query mass-to-charge ratio (m/z),
    find the peaks that are within the given PPM tolerance of the query m/z.
    """
    matching_peaks = []
    for peak in peaks:
        if mass_difference_in_ppm(mass1=peak.mz, mass2=query_mz) <= ppm_tolerance:
            matching_peaks.append(peak)
    return matching_peaks


def get_peak_product_ion_matches(
    spectrum: Spectrum,
    peptide: Union[Peptide, str],
    ion_types: Set[Literal[B_ION_TYPE, Y_ION_TYPE]] = {B_ION_TYPE, Y_ION_TYPE},
    peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
) -> List[PeakIonMatch]:
    """
    Compare the given spectrum to the given peptide. This method helps evaluate
    how strong the evidence is for a peptide-spectrum match (PSM).
    """
    if isinstance(peptide, str):
        peptide = Peptide(seq=peptide)
    # Product ions will have charge <= precursor's charge
    charges = list(range(1, spectrum.precursor_charge + 1))

    # Get product ions of the proposed peptide
    product_ions = peptide.product_ions(ion_types=ion_types, charges=charges)
    assert 2 * (len(peptide.seq) - 1) * len(charges) == len(product_ions)

    # Get peaks that match a product ion
    peak_ion_matches = []
    for ion in product_ions:
        matching_peaks = get_peaks_near_mz(
            query_mz=ion.mz,
            peaks=spectrum.peaks,
            ppm_tolerance=peak_to_ion_ppm_tolerance,
        )

        peak_ion_matches.extend(
            [
                PeakIonMatch(
                    ion_mz=ion.mz,
                    ion_charge=ion.charge,
                    ion_type=ion.ion_type,
                    ion_seq=ion.seq,
                    peak_mz=peak.mz,
                    peak_intensity=peak.intensity,
                    sample=spectrum.sample,
                    scan=spectrum.scan,
                )
                for peak in matching_peaks
            ]
        )

    return peak_ion_matches


def spectrum_peptide_plot(
    spectrum: Spectrum,
    seq: str,
    peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
) -> Axes:
    ion_intensity = max(peak.intensity for peak in spectrum.peaks) / 2
    peak_ion_matches = get_peak_product_ion_matches(
        spectrum=spectrum,
        peptide=seq,
        peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
    )
    peptide = Peptide(seq=seq)
    product_ions = peptide.product_ions(
        charges=list(range(1, spectrum.precursor_charge + 1)),
    )
    if ax is None:
        _, axs = fig_setup()
        ax = axs[0]

    # Plot spectrum
    spectrum.plot(
        ax=ax,
    )
    # Plot all product ions below spectrum
    plot_peaks(
        ax=ax,
        peaks=[Peak(mz=ion.mz, intensity=-1 * ion_intensity) for ion in product_ions],
    )
    # Plot ions that match peaks in a different color
    plot_peaks(
        ax=ax,
        peaks=[
            Peak(mz=match.ion_mz, intensity=-1 * ion_intensity)
            for match in peak_ion_matches
        ],
        color="red",
    )
    # Plot peaks that match ions in a different color
    plot_peaks(
        ax=ax,
        peaks=[
            Peak(mz=match.peak_mz, intensity=match.peak_intensity)
            for match in peak_ion_matches
        ],
        color="red",
        label="peak-to-ion matches",
    )
    # Plot y=0 line
    ax.axhline(0, color="black", linestyle="-", linewidth=0.5)

    # Finishing touches
    ax.set_ylim(bottom=-ion_intensity * 1.2)
    set_title_axes_labels(
        ax=ax,
        xlabel="m/z",
        ylabel="Intensity",
        title=(
            f"peptide: {seq}\nspectrum: {spectrum.uid}\nPPM tol: {peak_to_ion_ppm_tolerance}"
            if title is None
            else title
        ),
    )
    finalize(ax)
    return ax


@dataclass
class CometPSM:
    """Class for rows of Comet output"""

    sample: str
    num: int
    scan: int
    seq: str
    ions_matched: int
    ions_total: int
    proteins: List[str]
    precursor_charge: int
    spectrum_neutral_mass: float
    peptide_neutral_mass: float
    # retention_time: float
    # protein_count: int
    xcorr: float
    eval: float
    delta_cn: float
    q_value: Optional[float]

    @property
    def spectrum_uid(self) -> str:
        return Spectrum.get_uid(sample=self.sample, scan=self.scan)

    @property
    def uid(self) -> str:
        return self.spectrum_uid

    @property
    def is_hybrid(self):
        """
        A Comet PSM is a hybrid if the only proteins it appears in are hybrid proteins.
        If a PSM is in both a hybrid protein and a native protein, that means that the "hybrid"
        is a native sequence.
        """
        if all(self.check_if_hybrid_prot(prot=prot) for prot in self.proteins):
            return True
        else:
            return False

    @property
    def prop_ions_matched(self) -> float:
        return self.ions_matched / self.ions_total

    @staticmethod
    def check_if_hybrid_prot(prot: str):
        try:
            HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
            return True
        except ValueError:
            return False

    @classmethod
    def from_txt(
        cls,
        txt: str,
        as_df: bool = False,
        sample: str = "",
        only_num_1: bool = False,
        by_spectrum: bool = False,
    ) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads Comet results .txt file to a list of dataclasses or a dataframe
        """
        # Check whether TXT is from a direct Comet run or a Comet run via crux
        comet_txt = CometTxt(path=Path(txt))

        # Set sample if not provided
        if len(sample) == 0:
            sample = comet_txt.sample
        if comet_txt.file_type == CRUX:
            df = pd.read_csv(comet_txt.path, sep="\t")
            df[SAMPLE] = sample
            if "file" in df.columns:
                # If the 'file' column exists, it means it's the output of `crux assign-confidence`
                # in which case we need to set the sample differently
                df[SAMPLE] = df["file"].apply(
                    lambda file_path: Path(file_path).stem.split(".")[0]
                )
            # Convert crux column names to their Comet equivalent
            df.rename(
                columns={
                    "b/y ions matched": IONS_MATCHED,
                    "b/y ions total": IONS_TOTAL,
                    "xcorr score": XCORR,
                    "xcorr rank": NUM,
                    "protein id": PROTEIN,
                    "sequence": PLAIN_PEPTIDE,
                    "tdc q-value": Q_VAL,
                    "peptide mass": CALC_NEUTRAL_MASS,
                    "spectrum neutral mass": EXP_NEUTRAL_MASS,
                },
                inplace=True,
            )
        elif comet_txt.file_type == COMET:
            df = pd.read_csv(comet_txt.path, sep="\t", header=1)
            df[SAMPLE] = sample

        if as_df:
            return df
        else:
            data = [
                cls(
                    sample=row[SAMPLE],
                    scan=row[SCAN],
                    num=row[NUM],
                    ions_matched=row[IONS_MATCHED],
                    ions_total=row[IONS_TOTAL],
                    # protein_count=row[PROTEIN_COUNT],
                    proteins=row[PROTEIN].split(COMET_PROTEIN_SEPARATOR),
                    seq=row[PLAIN_PEPTIDE],
                    xcorr=row[XCORR],
                    eval=row[EVAL],
                    delta_cn=row[DELTA_CN],
                    q_value=row.get(Q_VAL, None),  # q-value is optional
                    precursor_charge=row[CHARGE],
                    # retention_time=row[RETENTION_TIME_STR],
                    spectrum_neutral_mass=row[EXP_NEUTRAL_MASS],
                    peptide_neutral_mass=row[CALC_NEUTRAL_MASS],
                )
                for _, row in df.iterrows()
            ]
            if only_num_1:
                data = [psm for psm in data if psm.num == 1]
            if by_spectrum:
                return organize_by_spectrum_uid(data=data)
            else:
                return data

    @classmethod
    def from_txts(
        cls, txts: List[Union[str, Path]], by_spectrum: bool = False
    ) -> List["CometPSM"]:
        spectrum_to_psms = defaultdict(list)
        for txt in txts:
            for psm in cls.from_txt(txt=txt):
                spectrum_to_psms[psm.spectrum_uid].append(psm)
        spectrum_to_psms = dict(spectrum_to_psms)
        if by_spectrum:
            return spectrum_to_psms
        else:
            return flatten_list_of_lists(list(spectrum_to_psms.values()))

    @staticmethod
    def get_top_psms(psms: List["CometPSM"]) -> List["CometPSM"]:
        """Get only `num=1` PSMs"""
        return [psm for psm in psms if psm.num == 1]

    def get_spectrum(self, spectra_dir: Union[str, Path] = Path("data")) -> Spectrum:
        matching_files = []
        for file in spectra_dir.rglob("*"):
            # file.name[-5:]
            if (
                file.is_file()
                and file.name[-5:] == ".mzML"
                and Mzml(path=file).name == self.sample
            ):
                matching_files.append(file)
        assert (
            len(matching_files) == 1
        ), f"Expected exactly one matching mzML file for PSM with sample {self.sample}, found {len(matching_files)}"
        return Mzml(path=matching_files[0]).get_spectrum(scan=self.scan)

    def to_psm(
        self,
        spectrum: Spectrum,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> "PSM":
        return PSM(
            spectrum=spectrum,
            seq=self.seq,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            xcorr=self.xcorr,
            q_value=self.q_value,
            prop_ions_matched=self.prop_ions_matched,
            positions=self.proteins,
        )

    @classmethod
    def load(cls, path: Union[str, Path], by_uid: bool = False):
        data = load_json(path=path)
        if isinstance(data, list):
            data = [cls(**psm_dict) for psm_dict in data]
        elif isinstance(data, dict):
            data = cls(**data)
        if by_uid:
            return organize_by_spectrum_uid(data=data)
        else:
            return data

    @staticmethod
    def save_psms_to_json(psms: List["CometPSM"], path: Path):
        to_json(
            data=[asdict(psm) for psm in psms],
            path=path,
        )

    def to_dict(self):
        data = asdict(self)
        data["uid"] = self.uid
        return data

    def save_to_json(self, path: Path):
        to_json(
            data=self.to_dict(),
            path=path,
        )

    @staticmethod
    def to_df(psms: List["CometPSM"]):
        return pd.DataFrame([psm.to_dict() for psm in psms])


@dataclass
class PSM:
    spectrum: Spectrum
    seq: str
    positions: List[str]
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    xcorr: Optional[float] = None
    q_value: Optional[float] = None
    comet_ions_matched: Optional[int] = None
    comet_ions_total: Optional[int] = None

    # Properties
    @cached_property
    def peak_ion_matches(self):
        return get_peak_product_ion_matches(
            spectrum=self.spectrum,
            peptide=self.seq,
            peak_to_ion_ppm_tolerance=self.peak_to_ion_ppm_tol,
        )

    @property
    def prop_comet_ions_matched(self):
        if (self.comet_ions_matched is None) or (self.comet_ions_total is None):
            return None
        else:
            return self.comet_ions_matched / self.comet_ions_total

    @property
    def df(self):
        return list_to_df(pydantic_list=self.peak_ion_matches)

    @property
    def num_b_ions_supported(self):
        if self.df is None:
            return 0
        return len(
            self.df[self.df.ion_type == B_ION_TYPE].groupby(
                by=["ion_charge", "ion_seq"]
            )
        )

    @property
    def uid(self):
        return self.spectrum.uid

    @property
    def num_ions_supported(self):
        return self.num_b_ions_supported + self.num_y_ions_supported

    @property
    def num_y_ions_supported(self):
        if self.df is None:
            return 0
        return len(
            self.df[self.df.ion_type == Y_ION_TYPE].groupby(
                by=["ion_charge", "ion_seq"]
            )
        )

    @property
    def prefixes_supported(self) -> List[str]:
        return list(self.sequences_supported(ion_type=B_ION_TYPE))

    @property
    def suffixes_supported(self) -> List[str]:
        return list(self.sequences_supported(ion_type=Y_ION_TYPE))

    @property
    def intensity_supported(self):
        return sum(
            [peak_ion_match.peak_intensity for peak_ion_match in self.peak_ion_matches]
        )

    @property
    def prop_intensity_supported(self):
        return self.intensity_supported / self.spectrum.total_intensity

    @property
    def prop_prefixes_supported(self):
        return len(self.prefixes_supported) / len(self.seq)

    @property
    def prop_suffixes_supported(self):
        return len(self.suffixes_supported) / len(self.seq)

    @property
    def mz_ppm_diff(self):
        seq_mz = compute_peptide_precursor_mz(seq=self.seq, charge=self.spectrum.z)
        return mass_difference_in_ppm(mass1=seq_mz, mass2=self.spectrum.mz)

    # Class methods
    @classmethod
    def from_spectrum_and_comet_psm(
        cls,
        spectrum: Spectrum,
        psm: CometPSM,
        ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> "PSM":
        return cls(
            spectrum=spectrum,
            seq=psm.seq,
            peak_to_ion_ppm_tol=ppm_tol,
            xcorr=psm.xcorr,
            q_value=psm.q_value,
            positions=psm.proteins,
            comet_ions_matched=psm.ions_matched,
            comet_ions_total=psm.ions_total,
        )

    @classmethod
    def from_comet_psms(
        cls,
        uid_to_spectrum: Dict[str, Spectrum],
        psms: List[CometPSM],
        ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> List["PSM"]:
        return [
            cls.from_spectrum_and_comet_psm(
                spectrum=uid_to_spectrum[psm.uid], psm=psm, ppm_tol=ppm_tol
            )
            for psm in psms
        ]

    # Instance methods
    def sequences_supported(
        self,
        ion_type: Literal[B_ION_TYPE, Y_ION_TYPE],
    ) -> Set[str]:
        if len(self.peak_ion_matches) == 0:
            return set()
        return set(self.df.loc[self.df.ion_type == ion_type, "ion_seq"].unique())

    def to_row(self) -> Dict:
        return {
            "uid": self.uid,
            "seq": self.seq,
            "xcorr": self.xcorr,
            "q_value": self.q_value,
            "prop_prefixes_supported": self.prop_prefixes_supported,
            "prop_suffixes_supported": self.prop_suffixes_supported,
            "prop_intensity_supported": self.prop_intensity_supported,
            "prefixes_supported": self.prefixes_supported,
            "suffixes_supported": self.suffixes_supported,
            "mz_ppm_diff": self.mz_ppm_diff,
            "positions": self.positions,
        }

    def plot(self):
        return spectrum_peptide_plot(
            spectrum=self.spectrum,
            seq=self.seq,
            peak_to_ion_ppm_tolerance=self.peak_to_ion_ppm_tol,
        )


def convert_comet_psms_to_psms(
    uid_to_spectrum_map: Dict[str, Spectrum],
    comet_psms: List[CometPSM],
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
) -> List[PSM]:
    psms = []
    for comet_psm in comet_psms:
        PSM.from_spectrum_and_comet_psm(
            spectrum=uid_to_spectrum_map[comet_psm.uid],
            psm=comet_psm,
            ppm_tol=peak_to_ion_ppm_tol,
        )
        psm = PSM(
            spectrum=uid_to_spectrum_map[comet_psm.spectrum_uid],
            seq=comet_psm.seq,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            xcorr=comet_psm.xcorr,
            q_value=comet_psm.q_value,
            prop_ions_matched=comet_psm.prop_ions_matched,
            positions=comet_psm.proteins,
        )
        psms.append(psm)
    return psms


class ProteinAbundance(BaseModel):
    protein_counts: Counter
    psms: Optional[List[CometPSM]] = None

    @property
    def df(self):
        return pd.DataFrame(
            {
                "protein": prot,
                "count": cnt,
            }
            for prot, cnt in self.protein_counts.items()
        )

    @classmethod
    def from_comet_txt(
        cls, txt: Union[str, Path], q_val_thresh: float = DEFAULT_Q_THRESHOLD
    ):
        psms = CometPSM.from_txt(txt=txt)
        return cls.from_comet_psms(psms=psms, q_threshold=q_val_thresh)

    @classmethod
    def from_comet_psms(
        cls, psms: List[CometPSM], q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> "ProteinAbundance":
        quality_psms = get_high_confidence_psms(
            psms=psms, score=Q_VAL, threshold=q_threshold
        )
        all_comet_proteins = flatten_list_of_lists(
            [psm.proteins for psm in quality_psms]
        )
        protein_counts = Counter(all_comet_proteins)
        return cls(protein_counts=protein_counts, psms=quality_psms)

    def get_top_n_protein_names(
        self, n: int, with_cnts: bool = False
    ) -> Union[Set[str], Dict[str, int]]:
        most_common_proteins = {
            prot: cnt for prot, cnt in self.protein_counts.most_common(n)
        }
        if with_cnts:
            return most_common_proteins
        else:
            return set(most_common_proteins.keys())

    def get_top_n_proteins(self, n: int, fasta: Path | str) -> List[Peptide]:
        return Fasta(path=fasta).get_proteins_by_name(
            names=self.get_top_n_protein_names(n=n)
        )

    def get_proteins_with_at_least_n_psms(
        self, n: int, fasta: Optional[Path | str] = None
    ) -> List[Union[str, Peptide]]:
        prot_names = [prot for prot, cnt in self.protein_counts.items() if cnt >= n]
        if fasta is not None:
            return Fasta(path=fasta).get_proteins_by_name(names=prot_names)
        else:
            return prot_names

    def relative_protein_abundances(
        self, fasta_path: Union[Path, str]
    ) -> "ProteinAbundance":
        fasta = Fasta(path=fasta_path)
        prot_name_to_leng = {prot.name: len(prot.seq) for prot in fasta.proteins}
        prot_cnts = defaultdict(int)
        for prot_name, cnt in self.protein_counts.items():
            prot_cnts[prot_name] = cnt / prot_name_to_leng[prot_name]
        return ProteinAbundance(protein_counts=Counter(prot_cnts))

    def plot_sorted_prot_cnts(
        self,
        top_n_prots: Optional[int] = None,
        ax: Optional[Axes] = None,
        title: Optional[str] = None,
        out_path: Optional[Union[str, Path]] = None,
    ) -> Axes:
        # Define data
        if len(self.protein_counts) == 0:
            logger.info("There are no PSMs to plot!")
            return
        items = sorted(self.protein_counts.items(), key=lambda x: x[1], reverse=True)
        if top_n_prots is not None:
            items = items[:top_n_prots]
        keys, values = zip(*items)

        # Plot
        if ax is None:
            fig, axs = fig_setup(h=8, w=10)
            ax = axs[0]
        ax.scatter(range(len(keys)), values)
        ax.set_xticks(range(len(keys)), keys, rotation=90, fontsize=8)

        # Add counts above points
        for i, v in enumerate(values):
            ax.annotate(
                f"{v}",
                (i, v),
                textcoords="offset points",
                xytext=(0, 8),
                ha="center",
                fontsize=9,
            )

        set_title_axes_labels(
            ax=ax,
            title=title,
            xlabel="Protein",
            ylabel="PSM counts",
        )
        finalize(ax)
        if out_path:
            save_fig(
                path=out_path,
            )
        return ax

    def get_ab(self, protein: str) -> int:
        return self.protein_counts[protein]

    def get_rel_ab(self, protein: str) -> float:
        max_count = max(self.protein_counts.values())
        return self.protein_counts[protein] / max_count

    def to_json(self, path: Union[str, Path]):
        data = dict(
            sorted(self.protein_counts.items(), key=lambda kv: kv[1], reverse=True)
        )
        to_json(data=data, path=path)

    def plot_counts_vs_prot_length(
        self,
        fasta: Union[str, Path, Fasta],
        title: Optional[str] = None,
        out_path: Optional[Union[str, Path]] = None,
    ):
        if isinstance(fasta, (str, Path)):
            fasta = Fasta(path=fasta)
        df = []
        for prot, cnt in self.protein_counts.items():
            df.append(
                {
                    "prot": prot,
                    "cnt": cnt,
                    "len": len(fasta.protein_name_to_seq_map[prot]),
                }
            )
        df = pd.DataFrame(df)
        fig, axs = fig_setup()
        ax = axs[0]
        sns.scatterplot(x=df.len, y=df.cnt, s=7, ax=ax, label=f"n={df.shape[0]}")
        set_title_axes_labels(
            ax=ax,
            title=title,
            xlabel="Protein length",
            ylabel="Number of PSMs from protein",
        )
        finalize(ax)
        if out_path:
            save_fig(
                path=out_path,
            )
        return df, ax


def fit_xcorr_to_qval_interpolator(
    psms: Union[pd.DataFrame, List[CometPSM]],
    # ax: Optional[Axes] = None,
) -> PchipInterpolator:
    # Get data
    if isinstance(psms, pd.DataFrame):
        df = psms.copy()
    else:
        df = pd.DataFrame(
            [(psm.xcorr, psm.q_value) for psm in psms],
            columns=[XCORR, Q_VAL],
        )
    xy = df.drop_duplicates(subset=XCORR)
    xy.sort_values(by=XCORR, inplace=True)
    x = xy[XCORR].to_numpy()
    y = xy[Q_VAL].to_numpy()
    interpolator = PchipInterpolator(x, y)
    return interpolator


def add_qvalue_interpolator_to_xcorr_plot(
    q_value_psms: Union[List[CometPSM], pd.DataFrame],
    ax: Axes,
    q_threshold: Optional[float] = None,
    ylabel: str = "Native log10(q-value)",
):
    q_interpolator = fit_xcorr_to_qval_interpolator(
        psms=q_value_psms,
    )
    ax_copy = ax.twinx()
    xmin, xmax = ax.get_xlim()
    x_new = np.linspace(xmin, xmax, 500)
    _ = ax_copy.plot(x_new, q_interpolator(x_new), "r--", label="Native q-value")
    if q_threshold is not None:
        _ = ax_copy.axhline(
            y=q_threshold, color="red", linestyle="--", label=f"q={q_threshold}"
        )
    ax_copy.set_yscale("log")  # set y-axis to log10 scale
    ax_copy.set_ylabel(ylabel, color="tab:red")
    ax_copy.tick_params(axis="y", labelcolor="tab:red")


def score_histogram(
    psms_by_type: Dict[str, List[Any]],
    score: str,
    ax: Optional[Axes] = None,
) -> Axes:
    if ax is None:
        _, axs = fig_setup()
        ax = axs[0]
    for key, psms in psms_by_type.items():
        if isinstance(psms, pd.DataFrame):
            data = psms[score]
        else:
            try:
                data = [getattr(psm, score) for psm in psms]
            except:
                data = psms
        _ = sns.kdeplot(
            data,
            ax=ax,
            label=f"{key} (n = {len(data)})",
        )
    return ax


def create_xcorr_dists_plot(
    psms_by_type: Dict[str, List[CometPSM]],
    title: Optional[str] = None,
    q_interpolating_psms: Optional[List[CometPSM]] = None,
    ax: Optional[Axes] = None,
) -> Axes:
    if ax is None:
        fig, axs = fig_setup()
        ax = axs[0]
    _ = score_histogram(psms_by_type=psms_by_type, score=XCORR, ax=ax)
    if q_interpolating_psms is not None:
        add_qvalue_interpolator_to_xcorr_plot(
            q_value_psms=q_interpolating_psms,
            ax=ax,
        )
    set_title_axes_labels(ax=ax, xlabel=XCORR, ylabel="Density", title=title)
    finalize(ax)
    return ax


@dataclass
class CometRunAnalysis:
    targets: List[CometPSM]
    decoys: List[CometPSM] = field(default_factory=list)
    assign_conf: List[CometPSM] = field(default_factory=list)
    interpolate: bool = True

    def __post_init__(self):
        if self.interpolate:
            logger.info("Interpolating q-values for all PSMs based on assign_conf PSMs")
            assert (
                len(self.assign_conf) > 0
            ), "Must have assign_conf psms to interpolate q-values"
            interpolator = fit_xcorr_to_qval_interpolator(psms=self.assign_conf)
            targets = []
            for psm in self.targets:
                psm.q_value = float(interpolator(psm.xcorr))
                targets.append(psm)
            self.targets = targets

            decoys = []
            for psm in self.decoys:
                psm.q_value = float(interpolator(psm.xcorr))
                decoys.append(psm)
            self.decoys = decoys

    @cached_property
    def top_targets(self) -> List[CometPSM]:
        return [psm for psm in self.targets if psm.num == 1]

    @cached_property
    def uid_to_top_target(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.top_targets}

    @cached_property
    def top_decoys(self) -> List[CometPSM]:
        return [psm for psm in self.decoys if psm.num == 1]

    @cached_property
    def uid_to_top_target(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.top_targets}

    @cached_property
    def uid_to_top_decoy(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.top_decoys}

    @cached_property
    def uid_to_assign_conf(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.assign_conf}

    @cached_property
    def target_df(self):
        return pd.DataFrame([psm.to_dict() for psm in self.targets])

    @cached_property
    def decoy_df(self):
        return pd.DataFrame([psm.to_dict() for psm in self.decoys])

    # @cached_property
    # def top_target_df(self):
    #     return pd.DataFrame([psm.to_dict() for psm in self.top_targets])

    # @cached_property
    # def top_decoy_df(self):
    #     return pd.DataFrame([psm.to_dict() for psm in self.top_targets])

    @property
    def top_psm_df(self):
        df = []
        for uid in set(self.uid_to_top_target.keys()).union(
            set(self.uid_to_top_decoy.keys())
        ):
            target_psm = self.uid_to_top_target.get(uid, None)
            decoy_psm = self.uid_to_top_decoy.get(uid, None)
            df.append(
                {
                    "uid": uid,
                    "target_xcorr": float(target_psm.xcorr) if target_psm else None,
                    "decoy_xcorr": float(decoy_psm.xcorr) if decoy_psm else None,
                    "target_seq": target_psm.seq if target_psm else None,
                    "decoy_seq": decoy_psm.seq if decoy_psm else None,
                }
            )
        return pd.DataFrame(df)

    def xcorr_target_vs_decoy_scatterplot(
        self, ax: Optional[Axes] = None, title: Optional[str] = None
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        sns.scatterplot(
            x=self.top_psm_df["decoy_xcorr"],
            y=self.top_psm_df["target_xcorr"],
            s=7,
            ax=ax,
            label=f"n={self.top_psm_df.shape[0]}",
        )
        plot_line(ax=ax, label="y=x")
        set_title_axes_labels(
            ax=ax,
            title=title,
            xlabel="Top decoy xcorr",
            ylabel="Top target xcorr",
        )
        finalize(ax)
        return ax

    def xcorr_target_vs_decoy_jointplot(
        self,
        title: Optional[str] = None,
    ):

        p = sns.jointplot(
            x=self.top_psm_df.target_seq.apply(len),
            y=self.top_psm_df["target_xcorr"] - self.top_psm_df["decoy_xcorr"],
            s=7,
            marginal_ticks=True,
            label=f"n={self.top_psm_df.shape[0]}",
        )
        p.set_axis_labels(
            xlabel="Top target peptide length",
            ylabel="top target xcorr - top decoy xcorr",
        )
        p.fig.suptitle(title)

        return p

    def xcorr_target_and_decoy_distributions(
        self, title: Optional[str] = None, ax: Optional[Axes] = None
    ) -> Axes:
        create_xcorr_dists_plot(
            psms_by_type={
                "Top targets": self.top_targets,
                "Top decoys": self.top_decoys,
            },
            q_interpolating_psms=(
                self.assign_conf if len(self.assign_conf) > 0 else None
            ),
            title=title,
            ax=ax,
        )

    def get_protein_abundance(
        self, q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> ProteinAbundance:
        if len(self.assign_conf) == 0:
            logger.warning(
                "Trying to get ProteinAbundance object when there are no assign-confidence PSMs. So skipping"
            )
            return
        return ProteinAbundance.from_comet_psms(
            psms=self.assign_conf, q_threshold=q_threshold
        )

    def basic_analysis(
        self,
        name: str,
        q_threshold: float = DEFAULT_Q_THRESHOLD,
        out_dir: Optional[Union[str, Path]] = None,
        fasta: Optional[Union[str, Path]] = None,
    ):
        create_xcorr_dists_plot(
            psms_by_type={
                "Top targets": self.top_targets,
                "Top decoys": self.top_decoys,
            },
            out_path=out_dir / "xcorr.png",
            title=name,
            q_interpolating_psms=(
                self.assign_conf if len(self.assign_conf) > 0 else None
            ),
        )
        if len(self.assign_conf) > 0:
            prot_ab = ProteinAbundance.from_comet_psms(
                psms=self.assign_conf, q_threshold=q_threshold
            )
            accepted_psms = [
                psm for psm in self.assign_conf if psm.q_value <= q_threshold
            ]
            title = f"{name}\nNumber of accepted PSMs (q<={q_threshold}): {len(accepted_psms)}"
            prot_ab.plot_sorted_prot_cnts(
                top_n_prots=100,
                title=title,
                out_path=out_dir
                / f"protein_abundance_q<={q_threshold}_topNprots100.png",
            )
            prot_ab.plot_sorted_prot_cnts(
                top_n_prots=30,
                title=title,
                out_path=out_dir
                / f"protein_abundance_q<={q_threshold}_topNprots30.png",
            )
            if fasta:
                prot_ab.plot_counts_vs_prot_length(
                    fasta=fasta,
                    title=title,
                    out_path=out_dir / f"prot_ab_vs_prot_len_q<={q_threshold}.png",
                )
            prot_ab.to_json(
                path=out_dir / f"protein_abundance_q<={q_threshold}.json",
            )


def get_high_confidence_psms(
    psms: List[CometPSM],
    score: Literal[Q_VAL] = Q_VAL,
    threshold: float = DEFAULT_Q_THRESHOLD,
) -> List[CometPSM]:
    """ """
    # Get only `num=1` PSMs
    psms = CometPSM.get_top_psms(psms=psms)
    return list(filter(lambda psm: getattr(psm, score) <= threshold, psms))


def get_comet_psm_to_spectrum_comparison_df(
    psms: List[CometPSM], spectrum_uid_to_spectrum: Dict[str, Spectrum]
):
    data = []
    for psm in psms:
        spectrum = spectrum_uid_to_spectrum[psm.spectrum_uid]
        data.append(
            (
                spectrum.precursor_mz,
                spectrum.charge,
                spectrum.precursor_intensity,
                spectrum.retention_time,
                psm.xcorr,
                psm.q_value,
                psm.ions_matched,
                psm.ions_total,
            )
        )
    return pd.DataFrame(
        data=data,
        columns=[
            "precursor_mz",
            "charge",
            "precursor_abundance",
            "retention_time",
            "xcorr",
            "q_value",
            "ions_matched",
            "ions_total",
        ],
    )
