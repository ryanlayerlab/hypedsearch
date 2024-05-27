import computational_pipeline.database_generator as database_generator
from preprocessing import preprocessing_utils
from main import get_spectra_files
from lookups.utils import ppm_to_da
from preprocessing.merge_search import modified_match_masses
import matplotlib.pyplot as plt

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25

prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database_generator.build_database(prot_path)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/NOD2_E3'
spectra_files = get_spectra_files(spectra_path)
spectra, boundaries = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

def calc_peaks_we_can_force_charge(input_spectrum, ppm_tol):
    conversion_counter = 0
    conditional_conversion_counter = 0
    masses = input_spectrum.mz_values
    converted_b_prec = masses[-2]
    converted_y_prec = masses[-1]
    b_tol = ppm_to_da(converted_b_prec, ppm_tol)
    y_tol = ppm_to_da(converted_y_prec, ppm_tol)
    for mass in masses[:-2]:
        if mass > converted_y_prec + y_tol:
            conversion_counter = conversion_counter + 1
        elif mass > converted_b_prec + b_tol:
            conditional_conversion_counter = conditional_conversion_counter + 1
    
    return conversion_counter, conditional_conversion_counter

def find_total_peaks_to_force_charge(spectra, ppm_tolerance):
    total_converted, total_conditional_converted = 0,0
    for spectrum in spectra:
        converted, conditional_converted = calc_peaks_we_can_force_charge(spectrum, ppm_tolerance)
        total_converted, total_conditional_converted = total_converted + converted, total_conditional_converted + conditional_converted
    print("We can force the charge for",total_converted ,"peaks and can conditionally force the charge for", total_conditional_converted, "peaks")
    
find_total_peaks_to_force_charge(spectra, ppm_tolerance)


def guess_hybrids(matched_masses_b, matched_masses_y):
    # Want to look at the peaks and see if we can guess whether this will be a hybrid.
    # In the case of a natural hit, we should have clusters which have hits of all sizes
    # In the case of hybrids, if the junction is in the middle somewhere, then we won't be able to find the big hits in the database
    
    # First step is to print the distribution of hits
    b_masses, y_masses, b_hit_nums, y_hit_nums = [],[], [], []
    for mass in matched_masses_b:
        b_masses.append(mass)
        b_hit_nums.append(len(matched_masses_b[mass]))
    for mass in matched_masses_y:
        y_masses.append(mass)
        y_hit_nums.append(len(matched_masses_y[mass]))
        
    return b_masses, b_hit_nums, y_masses, y_hit_nums
        
def mass_guess_hybrids(spectra, proteins, max_pep_len, ppm_tolerance, Natural):
    all_b_masses, all_b_hit_nums, all_y_hit_nums, all_y_masses = [],[],[],[]
    for spectrum in spectra:
        matched_masses_b, matched_masses_y = modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, False)
        b_masses, b_hit_nums, y_masses, y_hit_nums = guess_hybrids(matched_masses_b, matched_masses_y)
        all_b_masses, all_y_masses = all_b_masses + b_masses, all_y_masses + y_masses
        all_b_hit_nums, all_y_hit_nums = all_b_hit_nums + b_hit_nums, all_y_hit_nums + y_hit_nums
    
    fig1, ax1 = plt.subplots()
    ax1.scatter(all_b_masses, all_b_hit_nums, color='red', label="b hits")
    ax1.scatter(all_y_masses, all_y_hit_nums, color='blue', label="y hits")
    plt.title('Number of hits for a given size')
    plt.xlabel('Mass size')
    plt.ylabel('Number of hits')
    plt.legend()

    plt.savefig("Number_of_hits_per_size_for_"+Natural)
    
naturals = spectra[0:3] + spectra[12:]
mass_guess_hybrids(naturals, proteins, max_pep_len, ppm_tolerance, "natural")

hybrids = spectra[4:12]
mass_guess_hybrids(hybrids, proteins, max_pep_len, ppm_tolerance, "hybrid")
#later make sqlite queries that forces precursor masses to be matched with their associated ions
