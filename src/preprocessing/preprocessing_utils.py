from file_io import spectra
from utils import ppm_to_da
from constants import DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, PROTON_MASS, WATER_MASS
from gen_spectra import get_total_sum, convert_precursor_to_ion, get_precursor, gen_spectrum
from sqlite import database_file

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
    spectra_file, ppm_tol: int, peak_filter: int = 0, relative_abundance_filter: float = 0.0):
    linear_spectra = []
    all_spectra = []
    these_spectra = spectra.load(
        spectra_file, 
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    all_spectra += these_spectra
    # leave next 2 lines commented; uncomment only to test just specific indices
    # index_list = [3308] #For using a condensed database
    # these_spectra, all_spectra = reduce_database(all_spectra, these_spectra, index_list) 
    linear_spectra += list(set([
        x for spectrum in these_spectra for x in spectrum.mz_values
    ]))
    linear_spectra.sort()
    return (all_spectra)

def overlap_scoring(sequence, ppm_tol, input_masses):
    total_score = 0
    spectrum = gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tol)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical: #The bug is with 810 around here
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif abs(observed-theoretical) <= tol:
            total_score = total_score + 1
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return(total_score)

def distscoring(sequence, prec, prec_charge):
    combined_precursor = get_precursor(sequence, prec_charge)
    dist = abs(combined_precursor - prec)
    return 1/dist

def find_sequence(hit, protein_list):
    pid = hit[5]
    start, end = hit[1], hit[2]
    prot_seq = protein_list[pid][1]
    return prot_seq[start:end]  

def find_by_precursor(spectrum, prec_tol, protein_list, ppm_tol):
    dbf = database_file(10, False)
    converted_b, _ = convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    conv_prec_tol = ppm_to_da(converted_b, prec_tol)
    prec_hits = dbf.query_mass(converted_b, conv_prec_tol)[0]
    
    scored = []
    for hit in prec_hits:
        seq = find_sequence(hit, protein_list)
        distscore = distscoring(seq, spectrum.precursor_mass, spectrum.precursor_charge)
        aligned_score = overlap_scoring(seq, ppm_tol, spectrum.mz_values)
        scored.append((aligned_score, distscore, seq, hit))
        
    scored = sorted(scored, key = lambda x: (x[0], x[1]), reverse=True)
    return scored

def arePermutation(str1, str2):
    # Get lengths of both strings
    n1 = len(str1)
    n2 = len(str2)
 
    # If length of both strings is not same,
    # then they cannot be Permutation
    if (n1 != n2):
        return False
 
    # Sort both strings
    a = sorted(str1)
    str1 = " ".join(a)
    b = sorted(str2)
    str2 = " ".join(b)
 
    # Compare sorted strings
    for i in range(0, n1, 1):
        if (str1[i] != str2[i]):
            return False
 
    return True

def check_for_good_hit(scored_hits):
    top_hit = scored_hits[0]
    next_hit = scored_hits[1]
    counter = 2
    while arePermutation(top_hit[2], next_hit[2]):
        next_hit = scored_hits[counter]
        counter = counter + 1
    if top_hit[0] - next_hit[0] > 1: #the constant here can be changed with better data
        return True
    else:
        return False
