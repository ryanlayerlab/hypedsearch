from file_io import spectra
from utils import ppm_to_da

# turn the all spectra list into a list of boundaries
def make_boundaries(mz, ppm_tol):
    da_tol = ppm_to_da(mz, ppm_tol)
    return [mz - da_tol, mz + da_tol]

def load_spectra(
    spectra_files: list, ppm_tol: int, peak_filter: int = 0, relative_abundance_filter: float = 0.0):
    linear_spectra = []
    all_spectra = []
    for spectra_file in spectra_files:
        these_spectra = spectra.load(
            spectra_file, 
            peak_filter=peak_filter, 
            relative_abundance_filter=relative_abundance_filter
        )
        for object in these_spectra:
            object.mz_values.append(object.precursor_mass)
            object.abundance.append(100000000) #I gave it a crazy high abundance to represent precursor. Still a hack
        all_spectra += these_spectra
        # these_spectra = [these_spectra[0]]
        # all_spectra = [all_spectra[0]]
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.mz_values
        ]))
    linear_spectra.sort()
    boundaries = dict()
    for mz in linear_spectra:
        boundaries[mz] = make_boundaries(mz, ppm_tol)
    return (all_spectra, boundaries)