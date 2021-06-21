from file_io import spectra
from utils import ppm_to_da, make_overlap_boundaries
from dataclasses import dataclass, field

@dataclass
class Load_Spectra_Arguments:
    spectra_files: list = [], 
    ppm_tol: int = 0, 
    peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0

@dataclass
class Load_Spectra_Results:
    all_spectra: list = field(default_factory=list), 
    boundaries: list = field(default_factory=list), 
    mz_mapping: dict = field(default_factory=dict)

def load_spectra_file_into_memory(spectra_files: list = [], peak_filter: int = 0,relative_abundance_filter: float = 0.0):
    all_spectra = []
    linear_spectra = []
    for spectra_file in spectra_files:
        current_spectra = spectra.load(spectra_file, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        all_spectra += current_spectra
        linear_spectra += list(set([ x for spectrum in current_spectra for x in spectrum.spectrum]))
    return all_spectra, linear_spectra

def make_boundaries(mz,ppm_tol: int = 0):
    da_tol = ppm_to_da(mz, ppm_tol)
    return [mz - da_tol, mz + da_tol]

def make_mapping_for_mz_to_boundaries(linear_spectra,boundaries):
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(linear_spectra):
        if boundaries[b_i][0] <= linear_spectra[s_i] <= boundaries[b_i][1]:
            mz_mapping[linear_spectra[s_i]] = b_i 
            s_i += 1
        elif linear_spectra[s_i] < boundaries[b_i][0]:
            s_i += 1
        elif linear_spectra[s_i] > boundaries[b_i][1]:
            b_i += 1
    return mz_mapping

def load_spectra(lsa:Load_Spectra_Arguments)-> (Load_Spectra_Results):
    all_spectra,linear_spectra = load_spectra_file_into_memory(lsa.spectra_files,lsa.peak_filter,lsa.relative_abundance_filter)
    linear_spectra.sort()
    boundaries = [make_boundaries(mz,lsa.ppm_tol) for mz in linear_spectra]
    overlapped_boundaries = make_overlap_boundaries(boundaries)
    mz_mapping = make_mapping_for_mz_to_boundaries(linear_spectra,overlapped_boundaries)
    lsr = Load_Spectra_Results(all_spectra=all_spectra, boundaries=overlapped_boundaries, mz_mapping=mz_mapping)
    return lsr