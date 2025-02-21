from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from pydantic import BaseModel
from pyteomics import mzml


class Peak(BaseModel):
    mz: float
    abundance: float
    id: Optional[int] = None


class Spectrum(BaseModel):
    peaks: Tuple[Peak, ...]
    precursor_mz: float
    precursor_charge: int
    precursor_abundance: float
    spectrum_id: str
    retention_time: float
    mzml: Optional[str] = None
    scan_num: Optional[int] = None
    # For keeping track of whether the peaks have been processed:
    peaks_preprocessed: bool = False

    @classmethod
    def from_dict(cls, spectrum: Dict, mzml: Optional[str] = None) -> Spectrum:
        masses, abundances = tuple(spectrum["m/z array"]), tuple(
            spectrum["intensity array"]
        )
        peaks = [
            Peak(mz=masses[idx], abundance=abundances[idx], id=idx)
            for idx in range(len(masses))
        ]
        spectrum_id = spectrum.get("id", "")
        scan_num = int(spectrum_id.split("=")[1])
        return Spectrum(
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
    def from_mzml(cls, mzml_path: Union[str, Path]) -> List[Spectrum]:
        mzml_path = Path(mzml_path)
        spectra = mzml.read(str(mzml_path))
        return [
            cls.from_dict(spectrum=spectrum, mzml=mzml_path.stem)
            for spectrum in spectra
        ]

    def filter_to_top_n_peaks(self, n: int) -> None:
        # Update peaks
        new_peaks = top_n_peak_filtering(peaks=self.peaks, n=n)
        self.peaks = new_peaks

        # Update boolean that tracks whether peaks where preprocessed
        self.peaks_preprocessed = True


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
    abundances = [peak.abundance for peak in peaks]
    indices = get_indices_of_largest_elements(array=abundances, top_n=n)
    peaks = np.array(peaks)
    return list(peaks[indices])
