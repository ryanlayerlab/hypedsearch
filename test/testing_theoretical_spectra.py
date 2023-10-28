import gen_spectra
from utils import ppm_to_da
from scoring import scoring
from preprocessing.clustering import find_sequence
import database

def simulated_extensions(extended_b, extended_y, obs_prec, tol, precursor_charge, b_score, y_score):
    extensions = []
                    
    for b in extended_b:
        for y in extended_y:
            b_mass, y_mass = b[0], y[0]
            b_charge, y_charge = b[4], y[4]
            this_prec = gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, b_charge, y_charge, precursor_charge)
            if abs(this_prec - obs_prec) <= tol:
                return_tuple = (b + (b_score,),y + (y_score,))
                return return_tuple
    return extensions

def overlap_scoring(sequence, ppm_tol, input_masses):
    total_score = 0
    spectrum = gen_spectra.gen_spectrum(sequence)
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
    combined_precursor = gen_spectra.get_precursor(sequence, prec_charge)
    dist = abs(combined_precursor - prec)
    return 1/dist

def testing_score(input_masses, hypedsequence, specsequence, obs_prec, prec_charge, ppm_tol):
    hypedsearch_distance = distscoring(hypedsequence, obs_prec, prec_charge)
    hypedsearch_score = overlap_scoring(hypedsequence, ppm_tol, input_masses)
    spectrummill_distance = distscoring(specsequence, obs_prec, prec_charge)
    spectrummill_score = overlap_scoring(specsequence, ppm_tol, input_masses)
    if hypedsearch_score >= spectrummill_score:
        print("Hypedsearch scores just as well", hypedsearch_score, spectrummill_score)
    else:
        print("Hypedsearch scores worse than SpecMill", hypedsearch_score, spectrummill_score)
        
    if hypedsearch_distance > spectrummill_distance:
        print("Hypedsearch was closer to precursor", hypedsearch_distance, spectrummill_distance)
    else:
        print("SpectrumMill seq was closer to precursor", hypedsearch_distance, spectrummill_distance)    

def seq_score_by_dist(hybrid, prec_charge, b_pid, b_start, b_end, protein_list, y_pid, y_start, y_end):
    if hybrid:
        b_sequence = find_sequence(b_pid, b_start, b_end, protein_list)
        y_sequence = find_sequence(y_pid, y_start, y_end, protein_list)
        sequence = b_sequence + y_sequence
    else:
        sequence = find_sequence(b_pid, b_start, y_end, protein_list)
    combined_precursor = gen_spectra.get_precursor(sequence, prec_charge)
    dist = abs(combined_precursor - prec)
    return 1/dist

def testing_hypedsearch_score(extension, spectrum, ppm_tol, protein_list, max_len, prec, prec_charge, hybrid):
    b_pid, y_pid = extension[0][5], extension[1][5]
    b_start, b_end = extension[0][1], extension[0][2]
    y_start, y_end = extension[1][1], extension[1][2]
    if hybrid:
        distance_as_hybrid = scoring.score_by_dist(extension, prec, prec_charge, max_len, 1)
        score_as_hybrid = scoring.rescore(extension, spectrum, ppm_tol, 1, protein_list)
        score_as_hybrid_with_seq = seq_score_by_dist(hybrid, prec_charge, b_pid, b_start, b_end, protein_list, y_pid, y_start, y_end)
        print("Since this alignment is scored as a hybrid peptide, it has a score and distance score of:", score_as_hybrid, 1/distance_as_hybrid)
        print("If we score distance by finding seq first the distance score becomes", score_as_hybrid_with_seq)
    else:
        distance_as_natural = scoring.score_by_dist(extension, prec, prec_charge, max_len, 0)
        score_as_natural = scoring.rescore(extension, spectrum, ppm_tol, 0, protein_list)
        score_as_natural_with_seq = seq_score_by_dist(hybrid, prec_charge, b_pid, b_start, b_end, protein_list, y_pid, y_start, y_end)
        print("Since this alignment is scored as a natural peptide, it has a score and distance score of:", score_as_natural, 1/distance_as_natural)
        print("If we score distance by finding seq first the distance score becomes", score_as_natural_with_seq)

path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(path)

b_score = 4
y_score = 1
extended_b = [(286.157940435, 73, 78, 0, 2, 274), (321.676497435, 73, 79, 0, 2, 274), (378.218529435, 73, 80, 0, 2, 274), (442.739825935, 73, 81, 0, 2, 274)]
extended_y = [(124.06806228500001, 57, 59, 1, 2, 274), (188.58935878500003, 56, 59, 1, 2, 274), (266.63991428500003, 55, 59, 1, 2, 274), (344.690469785, 54, 59, 1, 2, 274), (388.206483785, 53, 59, 1, 2, 274), (452.253965285, 52, 59, 1, 2, 274)]
prec = 349.709089
prec_charge = 4
hybrid = False

# extension = simulated_extensions(extended_b, extended_y, prec, 10, prec_charge, b_score, y_score)

# target_b = (229.11782836914062, 75, 77, 1, 1, 126)
# target_y = (465.7381896972656, 19, 28, 1, 2, 235)
# target_extension = simulate_target_extension(target_b, target_y, prec, 10, prec_charge)
hypedsearch_result = "SLVLGLDLLAVEPG" #also try SLRL-LDLLAVEPG
spectrummill_result = "DLKIIWNKTKH"
sample_spectrum = [86.09611511230469,87.09937286376953,88.03799438476562,129.1030731201172,156.07711791992188,159.091064453125,201.12362670898438,261.06024169921875,278.1952819824219,314.1815490722656,314.6852722167969,338.4578552246094,349.112548828125,349.7074279785156,349.96002197265625,350.2093200683594,389.9046936035156,390.2397155761719,398.71099853515625,407.2212219238281,407.7216491699219,427.6009826660156,463.76177978515625,464.2640380859375,470.3000183105469,689.405619151242,698.410901501242]

# testing_hypedsearch_score(extension, sample_spectrum, 20, proteins.proteins, 23, prec, prec_charge, hybrid)
    
testing_score(sample_spectrum, hypedsearch_result, spectrummill_result, prec, prec_charge, 20)