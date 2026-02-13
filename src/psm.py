import logging
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union

import pandas as pd
from matplotlib.axes import Axes
from pydantic import BaseModel

from src.constants import (
    B_ION_TYPE,
    COMET,
    COMET_PROTEIN_SEPARATOR,
    CRUX,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_THRESHOLD,
    DELTA_CN,
    EVAL,
    HS_PREFIX,
    IONS_MATCHED,
    IONS_TOTAL,
    NUM,
    PLAIN_PEPTIDE,
    PROTEIN,
    Q_VAL,
    SAMPLE,
    SCAN,
    XCORR,
    Y_ION_TYPE,
)
from src.hybrids_via_clusters import HybridPeptide
from src.mass_spectra import Mzml, Peak, Spectrum, organize_by_spectrum_uid, plot_peaks
from src.peptides_and_ions import Fasta, Peptide, compute_peptide_precursor_mz
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
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
    path: Union[str, Path]
    sample: str = field(init=False)
    file_type: Literal["comet", "crux"] = field(init=False)
    # hybrid_run: Optional[bool]

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

    def read_psms(self, as_df: bool = False) -> Union[List["CometPSM"], pd.DataFrame]:
        """
        Reads the Comet output file to a list of dataclasses or a dataframe
        """
        if self.path is not None:
            return CometPSM.from_txt(txt=self.path, as_df=as_df)
        else:
            raise ValueError("No Comet output file found!")

    def get_header(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[1]

    def get_first_psm_line(self) -> str:
        if self.file_type == COMET:
            return self.path.read_text().split("\n")[2]

    def get_top_psms(self) -> List["CometPSM"]:
        psms = CometPSM.from_txt(txt=self.path)
        return [psm for psm in psms if psm.num == 1]


def read_comet_psms_from_dir(
    dir_path: Union[str, Path], glob_pattern: Optional[str] = None
):
    if glob_pattern is None:
        glob_pattern = "*.txt"
    all_psms = []
    for txt_file in Path(dir_path).glob(glob_pattern):
        psms = CometPSM.from_txt(txt=txt_file)
        all_psms.extend(psms)
    return all_psms


class PeakIonMatch(BaseModel):
    ion_mz: float
    ion_charge: int
    ion_seq: str
    ion_type: str
    peak_mz: float
    peak_intensity: float
    sample: Optional[str]
    scan: Optional[int]

    def mz_diff(self, type: Literal["rel", "rel_ppm"] = "rel_ppm"):
        """
        Mass (more precisely, m/z) difference between the theoretical ion and the peak.
        Let x_i = theoretical ion mass, x_p = peak mass,
        then returns
            - (x_i - x_t) / x_i when type='rel'
            - ((x_i - x_t) / x_i) * (10**6) when type='rel_ppm'
        """
        if type == "rel":
            return (self.ion_mz - self.peak_mz) / self.ion_mz
        elif type == "rel_ppm":
            return ((self.ion_mz - self.peak_mz) / self.ion_mz) * (10**6)

    @property
    def ion_id(self):
        return f"{self.ion_type}-{self.ion_seq}-z{self.ion_charge}"


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
    spectrum.plot_spectrum(
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
class PSM:
    spectrum: Spectrum
    seq: str
    positions: List[str]
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    xcorr: Optional[float] = None
    q_value: Optional[float] = None
    prop_ions_matched: Optional[float] = None

    @cached_property
    def peak_ion_matches(self):
        return get_peak_product_ion_matches(
            spectrum=self.spectrum,
            peptide=self.seq,
            peak_to_ion_ppm_tolerance=self.peak_to_ion_ppm_tol,
        )

    @property
    def df(self):
        return list_to_df(pydantic_list=self.peak_ion_matches)

    @property
    def num_b_ions_supported(self):
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
            df.rename(
                columns={
                    "b/y ions matched": IONS_MATCHED,
                    "b/y ions total": IONS_TOTAL,
                    "xcorr score": XCORR,
                    "xcorr rank": NUM,
                    "protein id": PROTEIN,
                    "sequence": PLAIN_PEPTIDE,
                    "tdc q-value": Q_VAL,
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
                    q_value=row.get(Q_VAL, None),  # Handle optional q-value
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
    ) -> PSM:
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
    def save_to_json(psms: List["CometPSM"], path: Path):
        to_json(
            data=[asdict(psm) for psm in psms],
            path=path,
        )

    def save(self, path: Path):
        to_json(
            data=asdict(self),
            path=path,
        )


def convert_comet_psms_to_psms(
    uid_to_spectrum_map: Dict[str, Spectrum],
    comet_psms: List[CometPSM],
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
) -> List[PSM]:
    psms = []
    for comet_psm in comet_psms:
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
        psms = get_high_confidence_psms(psms=psms, score=Q_VAL, threshold=q_threshold)
        all_comet_proteins = flatten_list_of_lists([psm.proteins for psm in psms])
        protein_counts = Counter(all_comet_proteins)
        return cls(protein_counts=protein_counts)

    def top_n_prots(
        self, n: int, with_cnts: bool = False
    ) -> Union[Set[str], Dict[str, int]]:
        most_common_proteins = {
            prot: cnt for prot, cnt in self.protein_counts.most_common(n)
        }
        if with_cnts:
            return most_common_proteins
        else:
            return set(most_common_proteins.keys())

    def relative_protein_abundances(
        self, fasta_path: Union[Path, str]
    ) -> "ProteinAbundance":
        fasta = Fasta(path=fasta_path)
        prot_name_to_leng = {prot.name: len(prot.seq) for prot in fasta.proteins}
        prot_cnts = defaultdict(int)
        for prot_name, cnt in self.protein_counts.items():
            prot_cnts[prot_name] = cnt / prot_name_to_leng[prot_name]
        return ProteinAbundance(protein_counts=Counter(prot_cnts))

    def plot(
        self, top_n_prots: Optional[int] = None, ax: Optional[Axes] = None
    ) -> Axes:
        # Define data
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
        set_title_axes_labels(
            ax=ax,
            # title="Protein counts",
            xlabel="Protein",
            ylabel="PSM counts",
        )
        finalize(ax)
        return ax

    def get_ab(self, protein: str) -> int:
        return self.protein_counts[protein]

    def get_rel_ab(self, protein: str) -> float:
        max_count = max(self.protein_counts.values())
        return self.protein_counts[protein] / max_count


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
