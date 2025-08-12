import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Dict, List, Optional, Union

import click
import numpy as np
from matplotlib.pyplot import Axes
from pydantic import BaseModel, BeforeValidator
from pyteomics import mzml as mzml_reader

from src.constants import SPECTRA_DIR, THOMAS_SAMPLES
from src.plot_utils import fig_setup, set_title_axes_labels
from src.utils import flatten_list_of_lists, to_path


@dataclass
class Peak:
    mz: float
    intensity: float
    id: Optional[int] = None


@dataclass
class Spectrum:
    precursor_mz: float
    precursor_charge: int
    precursor_abundance: float
    spectrum_id: str
    retention_time: float
    peaks: List[Peak] = field(default_factory=list, repr=False)
    mzml: Optional[Path] = None
    scan: Optional[int] = None
    # For keeping track of whether the peaks have been processed:
    peaks_preprocessed: bool = False

    @property
    def sample(self):
        if self.mzml is not None:
            return self.mzml.stem
        else:
            return None

    @property
    def total_intensity(self):
        return sum([peak.intensity for peak in self.peaks])

    @staticmethod
    def get_scan_number_from_id(spectrum_id: str) -> int:
        """
        Extract the scan number from the spectrum ID.
        """
        match = re.search(r"(?:scan|scanId)=(\d+)", spectrum_id)
        if match:
            return int(match.group(1))
        else:
            raise ValueError(f"Invalid spectrum ID format: {spectrum_id}")

    @classmethod
    def from_dict(cls, spectrum: Dict, mzml: Optional[Path] = None):
        # Extract scan number from 'id' key
        spectrum_id = spectrum.get("id")
        scan_num = int(re.search(r"(?:scan|scanId)=(\d+)", spectrum_id).group(1))

        # Get peaks
        masses, abundances = (
            tuple(spectrum["m/z array"]),
            tuple(spectrum["intensity array"]),
        )
        peaks = [
            Peak(mz=masses[idx], intensity=abundances[idx], id=idx)
            for idx in range(len(masses))
        ]
        return cls(
            scan=scan_num,
            peaks=peaks,
            spectrum_id=spectrum_id,
            mzml=mzml,
            precursor_mz=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
                "selectedIon"
            ][0]["selected ion m/z"],
            precursor_charge=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["charge state"],
            precursor_abundance=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["peak intensity"],
            retention_time=spectrum["scanList"]["scan"][0]["scan start time"],
        )

    @classmethod
    def parse_ms2_from_mzml(cls, spectra_file: Union[str, Path]) -> List["Spectrum"]:
        spectra_file = Path(spectra_file).absolute()
        ms2_spectra = []
        with mzml_reader.MzML(str(spectra_file)) as mzml:
            for spectrum in mzml:
                if spectrum["ms level"] == 1:
                    continue
                ms2_spectra.append(cls.from_dict(spectrum=spectrum, mzml=spectra_file))

            return ms2_spectra

    @classmethod
    def get_spectrum(cls, scan: int, mzml: Union[str, Path]):
        mzml = Path(mzml)
        with mzml_reader.MzML(str(mzml)) as reader:
            try:
                spectrum = reader.get_by_id(f"scan={scan}")
            except KeyError:
                spectrum = reader.get_by_id(f"scanId={scan}")
        return cls.from_dict(spectrum=spectrum, mzml=mzml)

    @classmethod
    def load_spectra_from_path(cls, path: Union[Path, str]) -> List["Spectrum"]:
        """
        Load spectra from a given path on the computer. If the path is a directory, find all
        mass spectra files in that directory and parse all their spectra
        """
        path = Path(path).absolute()
        if path.is_dir():
            spectra_files = list(path.glob("*.mzML"))
            spectra = [
                cls.parse_ms2_from_mzml(spectra_file=spectra_file)
                for spectra_file in spectra_files
            ]
            return flatten_list_of_lists(spectra)

        elif path.is_file():
            return cls.parse_ms2_from_mzml(spectra_file=path)

    def filter_to_top_n_peaks(self, n: int) -> None:
        if n > 0:
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
        color: str = "grey",
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        title = f"MZML={self.mzml.stem}; scan={self.scan}"
        plot_peaks(
            ax=ax,
            peaks=self.peaks,
            annotate=annotate,
            alpha=alpha,
            log_intensity=log_intensity,
            title=title,
            color=color,
        )
        ax.set_ylim(bottom=0)
        return ax


class Mzml(BaseModel):
    path: Annotated[Path, BeforeValidator(lambda x: to_path(path=x, check_exists=True))]

    @property
    def ms2_spectra(self) -> List["Spectrum"]:
        """
        Get all spectra from the mzML file.
        """
        return Spectrum.parse_ms2_from_mzml(spectra_file=self.path)

    @property
    def scans(self) -> List[int]:
        """
        Get all scan numbers from the mzML file.
        """
        return [spectrum.scan for spectrum in self.ms2_spectra]

    def get_spectrum(self, scan: int) -> "Spectrum":
        for spectrum in self.ms2_spectra:
            if spectrum.scan == scan:
                return spectrum

    @property
    def sample(self) -> str:
        return self.path.stem


def plot_peaks(
    ax: Axes,
    peaks: List[Peak],
    annotate: bool = True,
    log_intensity: bool = False,
    alpha: float = 1,
    color: str = "grey",
    label: str = "",
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
        spectra = Spectrum.parse_ms2_from_mzml(spectra_file=mzml_path)
        mzml_data.extend(list(spectra))
    return mzml_data


def get_spectrum_from_mzml(scan_num: int, mzml_path: Path):
    spectra = Spectrum.parse_ms2_from_mzml(spectra_file=mzml_path)

    spectrum = list(filter(lambda spectrum: spectrum.scan == scan_num, spectra))
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
    spectra = Spectrum.parse_ms2_from_mzml(spectra_file=mzml_path)
    matched_spectrum = None
    for spectrum in spectra:
        if spectrum.scan == scan_num:
            matched_spectrum = spectrum
            break
    return matched_spectrum


def get_mzml_for_sample(sample: str) -> Path:
    """
    Get the mzML file for a given sample.
    """
    mzml_path = SPECTRA_DIR / f"{sample}.mzML"
    if not mzml_path.exists():
        raise FileNotFoundError(f"MZML file for sample {sample} not found.")
    return mzml_path


def load_spectra_from_computer(path: Path) -> List[Spectrum]:
    """
    Load spectra from a given path on the computer. If the path is a directory, find all
    mass spectra files in that directory and parse all their spectra
    """
    if path.is_dir():
        spectra_files = list(path.glob("*.mzML"))
        spectra = [
            Spectrum.parse_ms2_from_mzml(spectra_file=spectra_file)
            for spectra_file in spectra_files
        ]
        return flatten_list_of_lists(spectra)

    elif path.is_file():
        return Spectrum.parse_ms2_from_mzml(spectra_file=path)

    spectra = Spectrum.parse_ms2_from_mzml(spectra_file=path)
    return spectra


@click.command(
    name="mzml-info",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help=("Print information about the given mzML file"),
)
@click.option(
    "--mzml",
    "-m",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the MZML file.",
)
def cli_mzml_info(mzml: Path):
    mzml = Mzml(path=mzml)
    print(f"MZML: {mzml.path}")
    print(f"\t - number of scans: {len(mzml.scans)}")
    if len(mzml.scans) < 50:
        print(f"\t - scans: {mzml.scans}")


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    cli.add_command(cli_mzml_info)
    cli()
