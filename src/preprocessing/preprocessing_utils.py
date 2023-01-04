from file_io import spectra
from utils import ppm_to_da
from constants import DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, PROTON_MASS, WATER_MASS
from gen_spectra import get_total_sum, convert_precursor_to_ion

# turn the all spectra list into a list of boundaries
def make_boundaries(mz, ppm_tol):
    da_tol = ppm_to_da(mz, ppm_tol)
    return [mz - da_tol, mz + da_tol]
def reduce_database(all_spectra, these_spectra, index_list):
    rall_spectra = []
    rthese_spectra = []
    for index in index_list:
        rall_spectra.append(all_spectra[index])
        rthese_spectra.append(these_spectra[index])
    return rall_spectra, rthese_spectra
    
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
            b_prec, y_prec = convert_precursor_to_ion(object.precursor_mass, object.precursor_charge)
            object.mz_values.append(b_prec)
            object.mz_values.append(y_prec)
            object.abundance.append(100000000) #I gave it a crazy high abundance to represent precursor. Still a hack
            object.abundance.append(100000000)
        all_spectra += these_spectra
        # leave commented; uncomment only to test just specific indices
        index_list = [5040] #For using a condensed database
        these_spectra, all_spectra = reduce_database(all_spectra, these_spectra, index_list) 
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.mz_values
        ]))
    linear_spectra.sort()
    boundaries = dict()
    for mz in linear_spectra:
        boundaries[mz] = make_boundaries(mz, ppm_tol)
    return (all_spectra, boundaries)