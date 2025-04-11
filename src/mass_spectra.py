from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from matplotlib.pyplot import Axes
from pydantic import BaseModel
from pyteomics import mzml

from src.constants import SPECTRA_DIR, THOMAS_SAMPLES
from src.plot_utils import fig_setup, set_title_axes_labels


@dataclass
class Peak:
    mz: float
    intensity: float
    id: Optional[int] = None


@dataclass
class Spectrum:
    peaks: Tuple[Peak, ...]
    precursor_mz: float
    precursor_charge: int
    precursor_abundance: float
    spectrum_id: str
    retention_time: float
    mzml: Optional[Path] = None
    scan_num: Optional[int] = None
    # For keeping track of whether the peaks have been processed:
    peaks_preprocessed: bool = False

    @classmethod
    def from_dict(cls, spectrum: Dict, mzml: Optional[Path] = None) -> "Spectrum":
        masses, abundances = tuple(spectrum["m/z array"]), tuple(
            spectrum["intensity array"]
        )
        peaks = [
            Peak(mz=masses[idx], intensity=abundances[idx], id=idx)
            for idx in range(len(masses))
        ]
        spectrum_id = spectrum.get("id", "")
        scan_num = int(spectrum_id.split("=")[1])
        return cls(
            # mass_over_charges=masses,
            # abundances=abundances,
            scan_num=scan_num,
            peaks=peaks,
            precursor_mz=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
                "selectedIon"
            ][0]["selected ion m/z"],
            precursor_charge=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["charge state"],
            precursor_abundance=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["peak intensity"],
            spectrum_id=spectrum_id,
            retention_time=spectrum["scanList"]["scan"][0]["scan start time"],
            mzml=mzml,
        )

    @classmethod
    def from_mzml(cls, mzml_path: Union[str, Path]) -> List["Spectrum"]:
        mzml_path = Path(mzml_path)
        spectra = mzml.read(str(mzml_path))
        return [
            cls.from_dict(spectrum=spectrum, mzml=mzml_path) for spectrum in spectra
        ]

    def filter_to_top_n_peaks(self, n: int) -> None:
        # Update peaks
        new_peaks = top_n_peak_filtering(peaks=self.peaks, n=n)
        self.peaks = new_peaks

        # Update boolean that tracks whether peaks where preprocessed
        self.peaks_preprocessed = True

    def plot_spectrum(
        self,
        ax: Optional[Axes] = None,
        annotate: bool = True,
        log_intensity: bool = False,
        alpha: float = 1,
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        title = f"MZML={self.mzml.stem}; scan={self.scan_num}"
        plot_peaks(
            ax=ax,
            peaks=self.peaks,
            annotate=annotate,
            alpha=alpha,
            log_intensity=log_intensity,
            title=title,
        )
        return ax


def plot_peaks(
    ax: Axes,
    peaks: List[Peak],
    annotate: bool = True,
    log_intensity: bool = False,
    alpha: float = 1,
    color: str = "grey",
    label: str = "peaks",
    title: Optional[str] = None,
):
    mzs = [peak.mz for peak in peaks]
    intensities = [peak.intensity for peak in peaks]
    if log_intensity:
        intensities = [np.log(intensity) for intensity in intensities]

    ax.vlines(
        mzs, [0], intensities, color=color, linewidth=0.5, alpha=alpha, label=label
    )
    if annotate:
        if log_intensity:
            ylabel = "log(intensity)"
        else:
            ylabel = "intensity"
        set_title_axes_labels(
            ax=ax,
            xlabel="m/z",
            ylabel=ylabel,
            title=title,
        )


def get_indices_of_largest_elements(array: List[float], top_n: int):
    if top_n >= len(array):
        return np.arange(0, len(array))
    array = np.array(array)
    # Get the indices of the largest N elements
    indices = np.argpartition(-array, top_n)[:top_n]
    # Sort these indices to have them in descending order of the elements
    sorted_indices = np.sort(indices)
    return sorted_indices


def top_n_peak_filtering(peaks: List[Peak], n: int) -> List[Peak]:
    abundances = [peak.intensity for peak in peaks]
    indices = get_indices_of_largest_elements(array=abundances, top_n=n)
    peaks = np.array(peaks)
    return list(peaks[indices])


def load_mzml_data(samples: List[str] = THOMAS_SAMPLES):
    mzml_data = []
    for sample in samples:
        print(f"Reading sample {sample}'s MZML")
        mzml_path = SPECTRA_DIR / f"{sample}.mzML"
        spectra = Spectrum.from_mzml(mzml_path=mzml_path)
        mzml_data.extend(list(spectra))
    return mzml_data


def get_spectrum_from_mzml(scan_num: int, mzml_path: Path):
    spectra = Spectrum.from_mzml(mzml_path=mzml_path)

    spectrum = list(filter(lambda spectrum: spectrum.scan_num == scan_num, spectra))
    assert (
        len(spectrum) == 1
    ), f"Scan number must be unique. There were {len(spectrum)} spectra with scan number {scan_num}."
    return spectrum[0]


def get_specific_spectrum_by_sample_and_scan_num(
    sample: Union[str, int], scan_num: int
) -> Spectrum:
    """
    Args:
        - scan_num is the spectrum's 1-based index
    """
    if isinstance(sample, str):
        mzml_path = SPECTRA_DIR / f"{sample}.mzML"
    elif isinstance(sample, int):
        mzml_path = SPECTRA_DIR / f"BMEM_AspN_Fxn{sample}.mzML"
    else:
        raise RuntimeError(
            f"Provided 'sample' should be type str or int. You provided {type(sample)}"
        )
    matched_spectrum = None
    spectra = Spectrum.from_mzml(mzml_path=mzml_path)
    matched_spectrum = None
    for spectrum in spectra:
        if spectrum.scan_num == scan_num:
            matched_spectrum = spectrum
            break
    return matched_spectrum
