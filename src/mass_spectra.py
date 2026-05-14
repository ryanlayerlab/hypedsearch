import logging
import re
import tempfile
from collections import Counter, defaultdict
from dataclasses import field
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import click
import numpy as np
import pandas as pd
import pymzml
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.pyplot import Axes
from pydantic import BaseModel
from pyteomics import mzml as mzml_reader

from src.constants import (
    COMMON_SPECTRA_ATTRS,
    DEFAULT_MAX_PRECURSOR_CHARGE,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    MAC_CRUX_EXECUTABLE,
    PRECURSOR_INTENSITY,
    SPECTRA_DF_NAME,
    SPECTRA_DIR,
    THOMAS_SAMPLES,
)
from src.plot_utils import (
    add_counts_to_histogram_boxes,
    fig_setup,
    finalize,
    plot_histogram,
    save_fig,
    set_title_axes_labels,
)
from src.utils import (
    CmdLineRunner,
    compute_gini_coefficient,
    flatten_list_of_lists,
    load_json,
    mass_difference_in_ppm,
    save_pydantic_objects_to_json,
)

logger = logging.getLogger(__name__)


class Peak(BaseModel):
    mz: float
    intensity: float
    id: Optional[int] = None


class Spectrum(BaseModel):
    precursor_mz: float
    precursor_charge: int
    precursor_intensity: float
    spectrum_id: str
    retention_time: float
    peaks: List[Peak] = field(default_factory=list, repr=False)
    mzml: Optional[Path] = None
    scan: Optional[int] = None
    mzml_index: Optional[int] = None
    # For keeping track of whether the peaks have been processed:
    peaks_preprocessed: bool = False

    @property
    def sample(self):
        if self.mzml is not None:
            return self.mzml.stem
        else:
            return None

    @property
    def uid(self):
        if self.sample is None:
            raise RuntimeError(
                f"Spectrum scan={self.scan} has no sample associated to it"
            )
        return self.get_uid(sample=self.sample, scan=self.scan)

    def get_precursor_peaks(
        self, peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        precursor_peaks = list(
            filter(
                lambda peak: abs(
                    mass_difference_in_ppm(mass1=peak.mz, mass2=self.precursor_mz)
                )
                <= peak_to_ion_ppm_tol,
                self.peaks,
            )
        )
        return precursor_peaks

    def get_non_precursor_peaks(
        self,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ):
        non_precursor_peaks = list(
            filter(
                lambda peak: not (
                    abs(mass_difference_in_ppm(mass1=peak.mz, mass2=self.precursor_mz))
                    <= peak_to_ion_ppm_tol
                ),
                self.peaks,
            )
        )
        return non_precursor_peaks

    def get_total_intensity(
        self, peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ) -> float:
        return sum(
            [
                peak.intensity
                for peak in self.get_non_precursor_peaks(
                    peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
            ]
        )

    @property
    def mz(self):
        return self.precursor_mz

    @property
    def z(self):
        return self.precursor_charge

    @property
    def charge(self):
        return self.precursor_charge

    @property
    def noise(self):
        pass

    @staticmethod
    def get_uid(sample: str, scan: int) -> str:
        return f"mzml={sample};scan={scan}"

    @staticmethod
    def parse_uid(uid: str) -> Tuple[str, int]:
        match = re.match(r"^mzml=(.+);scan=(\d+)$", uid)
        if match:
            sample = match.group(1)
            scan = int(match.group(2))
            return sample, scan
        else:
            raise ValueError(f"Invalid UID format: {uid}")

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
        scan_num = cls.get_scan_number_from_id(spectrum_id=spectrum_id)

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
            mzml_index=spectrum.get("index"),
            peaks=peaks,
            spectrum_id=spectrum_id,
            mzml=mzml,
            precursor_mz=spectrum["precursorList"]["precursor"][0]["selectedIonList"][
                "selectedIon"
            ][0]["selected ion m/z"],
            precursor_charge=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["charge state"],
            precursor_intensity=spectrum["precursorList"]["precursor"][0][
                "selectedIonList"
            ]["selectedIon"][0]["peak intensity"],
            retention_time=spectrum["scanList"]["scan"][0]["scan start time"],
        )

    @classmethod
    def parse_ms2_from_mzml(
        cls, mzml: Union[str, Path], by_uid: bool = False, as_df: bool = False
    ) -> Union[List["Spectrum"], Dict[str, "Spectrum"], pd.DataFrame]:
        logger.info(f"Parsing MS2 spectra from mzML file: {mzml}")
        mzml_path = Path(mzml).absolute()
        ms2_spectra = []
        with mzml_reader.MzML(str(mzml_path)) as mzml:
            for spectrum in mzml:
                if spectrum["ms level"] != 2:
                    continue
                spectrum = cls.from_dict(spectrum=spectrum, mzml=mzml_path)
                ms2_spectra.append(spectrum)
        if by_uid:
            return organize_by_spectrum_uid(data=ms2_spectra)
        else:
            if as_df:
                return cls.to_df(spectra=ms2_spectra)
            else:
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
        `*.mzML` files in that directory and parses their spectra
        """
        path = Path(path).absolute()
        if path.is_dir():
            spectra_files = list(path.glob("*.mzML"))
            spectra = [
                cls.parse_ms2_from_mzml(mzml=spectra_file)
                for spectra_file in spectra_files
            ]
            return flatten_list_of_lists(spectra)

        elif path.is_file():
            return cls.parse_ms2_from_mzml(mzml=path)

    @staticmethod
    def to_df(spectra: List["Spectrum"]) -> pd.DataFrame:
        return create_spectra_df(spectra=spectra)

    def filter_to_top_n_peaks(self, n: int) -> None:
        if n > 0:
            # Update peaks
            new_peaks = top_n_peak_filtering(peaks=self.peaks, n=n)
            self.peaks = new_peaks

            # Update boolean that tracks whether peaks where preprocessed
            self.peaks_preprocessed = True

    @property
    def plot_title(self):
        return (
            f"MZML={self.mzml.stem}; scan={self.scan}\n"
            + f"RT={self.retention_time:.2f}; m/z={self.precursor_mz:.2f}; z={self.precursor_charge}; intensity={self.precursor_intensity:.1f}"
        )

    def plot(
        self,
        ax: Optional[Axes] = None,
        annotate: bool = True,
        log_intensity: bool = False,
        alpha: float = 1,
        color: str = "grey",
        peak_to_ion_ppm_tol: Optional[float] = None,
    ):
        """
        If peak_to_ion_ppm_tol is given, indicate precursor peak with a different color.
        """
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        plot_peaks(
            ax=ax,
            peaks=self.peaks,
            annotate=annotate,
            alpha=alpha,
            log_intensity=log_intensity,
            title=self.plot_title,
            color=color,
        )
        top = None
        if peak_to_ion_ppm_tol is not None:
            precursor_peaks = self.get_precursor_peaks(
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
            )
            if len(precursor_peaks) > 0:
                plot_peaks(
                    ax=ax,
                    peaks=precursor_peaks,
                    color="black",
                    label="precursor peak",
                    lw=1,
                )
                # top = (
                #     max(
                #         peak.intensity
                #         for peak in self.get_non_precursor_peaks(
                #             peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                #         )
                #     )
                #     * 1.1
                # )
        finalize(ax)
        ax.set_ylim(bottom=0, top=top)
        return ax

    @staticmethod
    def plot_charges(
        spectra: List["Spectrum"],
        ax: Optional[Axes] = None,
        title: Optional[str] = None,
        out_path: Optional[Union[str, Path]] = None,
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        _ = sns.histplot([sp.precursor_charge for sp in spectra], ax=ax)
        add_counts_to_histogram_boxes(ax=ax)
        set_title_axes_labels(
            ax=ax, title=title, xlabel="precursor charge", ylabel="Count"
        )
        finalize(axs)
        if out_path is not None:
            save_fig(out_path)

    @staticmethod
    def plot_spectra_info(
        spectra: List["Spectrum"],
        attrs: List[str] = COMMON_SPECTRA_ATTRS,
        add_counts: bool = True,
    ) -> Tuple[Figure, List[Axes]]:
        fig, axs = fig_setup(nrows=len(attrs), ncols=1)
        for i, attr in enumerate(attrs):
            if attr == "precursor_intensity":
                plot_histogram(
                    objects=spectra,
                    attr=attr,
                    ax=axs[i],
                    add_counts=add_counts,
                    rug_plot=True,
                )
            else:
                plot_histogram(
                    objects=spectra, attr=attr, ax=axs[i], add_counts=add_counts
                )
        finalize(axs)
        return fig, axs

    @classmethod
    def load_spectra_from_mzmls(cls, mzmls: List[Union[str, Path]]) -> List["Spectrum"]:
        all_spectra = []
        for mzml in mzmls:
            spectra = cls.parse_ms2_from_mzml(mzml=mzml)
            all_spectra.extend(spectra)
        return all_spectra

    @staticmethod
    def save_to_json(spectra: List["Spectrum"], path: Path):
        save_pydantic_objects_to_json(
            objects=spectra,
            path=path,
        )

    def to_dict(self):
        return self.model_dump(mode="json")

    @classmethod
    def load_from_json(
        cls, json: Union[str, Path], by_uid: bool = False
    ) -> List["Spectrum"]:
        data = load_json(path=json)
        spectra = [cls(**d) for d in data]  # took 14s
        if by_uid:
            return organize_by_spectrum_uid(data=spectra)
        else:
            return spectra

    def number_of_peaks_needed_to_capture_percent_of_total_intensity(
        self, percent: float
    ) -> int:
        """
        Get the number of peaks needed to capture the given percent of total intensity.
        """
        prop_of_tot_intensity = (
            np.cumsum(sorted([peak.intensity for peak in self.peaks], reverse=True))
            / self.total_intensity
        )
        num_peaks = np.searchsorted(prop_of_tot_intensity, percent, side="right") + 1
        return num_peaks

    def get_gini_coefficient(self, peak_to_ion_ppm_tol: float) -> float:
        return compute_gini_coefficient(
            values=[
                peak.intensity
                for peak in self.get_non_precursor_peaks(
                    peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
            ]
        )

    @property
    def info(
        self,
        #  percent: float = 0.75
    ) -> Dict:
        return {
            "uid": self.uid,
            "mz": self.precursor_mz,
            "z": self.precursor_charge,
            "intensity": self.precursor_intensity,
            "rt": self.retention_time,
            "num_peaks": len(self.peaks),
            # "total_intensity": self.total_intensity,
            "sample": self.sample,
            "scan": self.scan,
            # "gini": self.gini,
            # f"num_peaks_to_capture{percent}": self.number_of_peaks_needed_to_capture_percent_of_total_intensity(
            # percent=percent
            # ),
        }

    @staticmethod
    def create_spectra_df(spectra: List["Spectrum"]):
        df = pd.DataFrame(
            data=[spectrum.info for spectrum in spectra],
        )
        return df


def organize_by_spectrum_uid(data: List[Any]):
    if len(data) == 0:
        return {}
    else:
        spec_id = getattr(data[0], "spectrum_uid", None)
        if spec_id is not None:
            uid_to_objs = defaultdict(list)
            for datum in data:
                uid_to_objs[datum.spectrum_uid].append(datum)
            return dict(uid_to_objs)
        else:
            # So it can work on spectra too
            return {datum.uid: datum for datum in data}


class Mzml(BaseModel):
    path: Path

    @cached_property
    def ms2_spectra(self) -> List["Spectrum"]:
        """
        Get all spectra from the mzML file.
        """
        return Spectrum.parse_ms2_from_mzml(mzml=self.path)

    @cached_property
    def id_to_spectrum(self) -> Dict[str, Spectrum]:
        return {spectrum.uid: spectrum for spectrum in self.ms2_spectra}

    @property
    def scan_numbers(self) -> List[int]:
        """
        Get all scan numbers from the mzML file.
        """
        scan_numbers = []
        with pymzml.run.Reader(str(self.path)) as reader:
            for spec in reader:
                if spec.ms_level > 1:  # MS2 or higher
                    scan_numbers.append(int(spec.ID))
        return scan_numbers

    def get_spectrum(self, scan: int) -> "Spectrum":
        return Spectrum.get_spectrum(scan=scan, mzml=self.path)

    @property
    def sample(self) -> str:
        return self.path.stem

    @property
    def name(self) -> str:
        return self.get_mzml_name(mzml=self.path)

    @staticmethod
    def get_mzml_name(mzml: Union[str, Path]):
        """Remove .mzML extension"""
        mzml = Path(mzml)
        return f"{mzml.name[:-5]}"

    def run_param_medic(
        self,
        out_dir: Optional[str | Path] = None,
        crux_path: str | Path = MAC_CRUX_EXECUTABLE,
    ) -> pd.DataFrame:
        param_medic_suffix = "param-medic.txt"
        create_cmd = lambda output_dir: [
            f"{crux_path} param-medic",
            f'"{self.path}"',
            "--overwrite T",
            f'--fileroot "{self.name}"',
            f'--output-dir "{output_dir}"',
        ]
        if out_dir:
            result = CmdLineRunner.run_cmd(cmd=create_cmd(output_dir=out_dir))
            df = pd.read_csv(
                Path(out_dir) / f"{self.name}.{param_medic_suffix}", sep="\t"
            )
        else:
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_dir = Path(tmp_dir)
                result = CmdLineRunner.run_cmd(cmd=create_cmd(output_dir=tmp_dir))
                df = pd.read_csv(
                    Path(tmp_dir) / f"{self.name}.{param_medic_suffix}", sep="\t"
                )
        return result, df


def plot_peaks(
    ax: Axes,
    peaks: List[Peak],
    annotate: bool = True,
    log_intensity: bool = False,
    lw: float = 0.5,
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
        mzs, [0], intensities, color=color, linewidth=lw, alpha=alpha, label=label
    )
    # Plot y=0 line
    ax.axhline(0, color="black", linestyle="-", linewidth=0.5)

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


def spectra_pairplot(spectra: List[Spectrum]):
    spec_df = pd.DataFrame(
        data=[
            (sp.precursor_mz, sp.charge, sp.precursor_intensity, sp.retention_time)
            for sp in spectra
        ],
        columns=["precursor_mz", "charge", "precursor_abundance", "retention_time"],
    )
    g = sns.pairplot(
        data=spec_df,
        corner=True,
        plot_kws={"s": 7},
        diag_kind="hist",
    )
    for i in range(len(g.axes)):
        g.axes[i, i].set_visible(False)


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
        spectra = Spectrum.parse_ms2_from_mzml(mzml=mzml_path)
        mzml_data.extend(list(spectra))
    return mzml_data


def get_mzml_for_sample(sample: str) -> Path:
    """
    Get the mzML file for a given sample.
    """
    mzml_path = SPECTRA_DIR / f"{sample}.mzML"
    if not mzml_path.exists():
        raise FileNotFoundError(f"MZML file for sample {sample} not found.")
    return mzml_path


def load_spectra_from_path(path: Union[str, Path]) -> List[Spectrum]:
    """
    Load spectra from a given path on the computer. If the path is a directory, find all
    mass spectra files in that directory and parse all their spectra
    """
    path = Path(path)
    if path.is_dir():
        spectra_files = list(path.glob("*.mzML"))
        spectra = [
            Spectrum.parse_ms2_from_mzml(mzml=spectra_file)
            for spectra_file in spectra_files
        ]
        return flatten_list_of_lists(spectra)

    elif path.is_file():
        return Spectrum.parse_ms2_from_mzml(mzml=path)

    else:
        raise RuntimeError(f"Path {path} is neither a file nor a directory.")


def precursor_mz_plot(
    spectra: List[Spectrum],
    title: str,
) -> Figure:
    # Create data to plot
    mz_to_spectra = defaultdict(list)
    for spectrum in spectra:
        mz_to_spectra[spectrum.precursor_mz].append(spectrum)
    mz_to_spectra = dict(mz_to_spectra)
    df = pd.DataFrame(
        {
            "precursor_mz": mz,
            "num_spectra": len(spectra),
        }
        for mz, spectra in mz_to_spectra.items()
    )

    # Plot
    fig, axs = fig_setup()
    _ = sns.scatterplot(data=df, x="precursor_mz", y="num_spectra", ax=axs[0], s=7)
    set_title_axes_labels(
        ax=axs[0],
        xlabel="precursor m/z",
        ylabel="Number of spectra",
    )
    finalize(axs)
    # fig = interactive_scatter_plot(
    #     df=df, x_colm="precursor_mz", y_colm="num_spectra", title=title
    # )
    return fig


def plot_spectra_histograms(spectra: List[Spectrum], add_cnts: bool = False):
    df = create_spectra_df(spectra=spectra)
    attrs = COMMON_SPECTRA_ATTRS + ["gini"]
    fig, axs = fig_setup(nrows=len(attrs), ncols=1)
    for i, attr in enumerate(attrs):
        if attr == PRECURSOR_INTENSITY:
            # Include rug plot for precursor_intensity
            plot_histogram(
                values=df[attr], ax=axs[i], rug_plot=True, add_counts=add_cnts
            )
        else:
            plot_histogram(values=df[attr], ax=axs[i], add_counts=add_cnts)
    finalize(axs)
    return df, fig, axs


def create_spectra_df(
    spectra: List[Spectrum],
):
    rows = []
    for spectrum in spectra:
        data = {"uid": spectrum.uid}
        for attr in COMMON_SPECTRA_ATTRS:
            data[attr] = getattr(spectrum, attr)
            data["gini"] = compute_gini_coefficient(
                values=[peak.intensity for peak in spectrum.peaks]
            )
        rows.append(data)
    df = pd.DataFrame(rows)
    return df


def peak_intensity_captured_plot(spectra: List[Spectrum]):
    fig, axs = fig_setup(ncols=2)
    for spectrum in spectra:
        prop_of_tot_intensity = (
            np.cumsum(sorted([peak.intensity for peak in spectrum.peaks], reverse=True))
            / spectrum.total_intensity
        )
        x = np.array(range(1, len(prop_of_tot_intensity) + 1))
        _ = axs[0].plot(
            x,
            prop_of_tot_intensity,
            # "k-"
        )
        x = np.array(range(1, len(prop_of_tot_intensity) + 1)) / len(
            prop_of_tot_intensity
        )
        _ = axs[1].plot(
            x,
            prop_of_tot_intensity,
            # "k-"
        )
    set_title_axes_labels(
        ax=axs[0],
        xlabel="Number of peaks",
        ylabel="Proportion of total\nintensity captured",
    )
    set_title_axes_labels(
        ax=axs[1],
        xlabel="Proportion of peaks",
        ylabel="Proportion of total\nintensity captured",
    )
    finalize(axs)
    return fig, axs


def create_sample_scan_to_spectrum_map(
    spectra_dir: Optional[Union[str, Path]] = None, mzmls: Optional[List[Path]] = None
) -> Dict[Tuple[str, int], Spectrum]:
    sample_scan_to_spectrum_map = {}
    if mzmls is not None:
        for mzml in mzmls:
            for spectrum in Mzml(path=mzml).ms2_spectra:
                sample_scan_to_spectrum_map[spectrum.uid] = spectrum
    else:
        for mzml in Path(spectra_dir).glob("*.mzML"):
            for spectrum in Mzml(path=mzml).ms2_spectra:
                sample_scan_to_spectrum_map[spectrum.uid] = spectrum
    return sample_scan_to_spectrum_map


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
