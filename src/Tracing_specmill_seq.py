import database
import os
from multiset import *
import collections
import numpy as np
from preprocessing import preprocessing_utils, clustering, merge_search
from scoring.scoring import modified_losing_water, modified_losing_ammonium
import gen_spectra
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import modified_match_masses
from sqlite import database_file
from gen_spectra import get_raw_mass
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25

spectra_path = '/home/naco3124/snakemake/spectra/'
spectra_files = get_spectra_files(spectra_path)
spectra = []

def overlap_scoring(sequence, input_masses, ppm_tolerance):
    total_score = 0
    spectrum = gen_spectra.gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tolerance)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical:
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


def get_natives_and_hybrids(filepath):
    natives, hybrids = dict(), dict()
    with open(filepath, 'r') as truth_set:
        for i, line in enumerate(truth_set):
            split_line = line.split('\t')
            if split_line[21] == "DGLNHL":
                print("here")
            if "Hybrid" not in line:
                if split_line[21] not in natives.keys():
                    natives[split_line[21]] = []
                natives[split_line[21]].append(split_line[17])

            else:
                if split_line[21] not in hybrids.keys():
                    hybrids[split_line[21]] = []
                hybrids[split_line[21]].append(split_line[17])
            
    return natives.keys(), natives, hybrids.keys(), hybrids

def get_hypedsearch_results(file):
    result=[]
    with open(file, 'r') as f:
        prev_id='spectrum_id'
        prev_score, prev_abundance = -1, -1
        for line in f:
            A = line.split("\t")
            spec_id = A[0]
            score = A[3]
            abundance = A[4]
            if score == prev_score and abundance == prev_abundance:
                prev_score = score
                prev_abundance = abundance
                seq = A[2]
                eval = A[1]
                result.append((spec_id, score, seq, eval))
            elif spec_id != prev_id:
                prev_id = spec_id
                prev_score = score
                prev_abundance = abundance
                seq = A[2]
                eval = A[1]
                result.append((spec_id, score, seq, eval))
    return result

def get_target_score(input_masses, target, ppm_tol):
    target_score = overlap_scoring(target, input_masses, ppm_tol)
    minus_water, _, _ = modified_losing_water(target, input_masses, ppm_tol)
    minus_ammonium, _, _ = modified_losing_ammonium(target, input_masses, ppm_tolerance)
    return (target_score + minus_water + minus_ammonium) / len(target)

def get_output_files(spectra_folder):
    spectra_files = []
    output_dict = dict()
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    
    for file in spectra_files:
        filesplit = file.split("Fxn")[1]
        file_id = int(filesplit.split(".")[0])
        output_dict[file_id] = file
    return output_dict

output_dict = get_output_files("/home/naco3124/snakemake/output/Hypedsearch_outputs")

# count = 0;
# for file in spectra_files:
#     spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
#     filesplit = file.split("Fxn")[1]
#     file_id = filesplit.split(".")[0]
#     hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
#     for answer in hypedsearch_answers:
#         score = answer[1]
#         if int(score) == 0:
#             count = count + 1

# print("Hypedsearch returned", count, "seqs with score of 0")

with open("Sequence_tracing.txt", "w") as w:
    w.write("")
        
with open("Sequence_tracing.txt", "a") as w:
    
    w.write("When Hypedsearch missed a target hybrid:\n\n")
    
    for file in spectra_files:
        filename = file.split("/")[-1]
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]

        hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
        for result in hypedsearch_answers:
            spec_id = result[0]
            hypedsearch_score = result[1]
            seq = result[2]
            hybrid_eval = result[3]
            input_spectrum = spectra[int(spec_id)]
            prec_tolerance = ppm_to_da(input_spectrum.precursor_mass, prec_tol)
            
            if abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALEVE", ppm_tolerance)
                if tscore > float(hypedsearch_score):
                    w.write("For " + filename + " DLQTLAL-EVE scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALNAAR", ppm_tolerance)
                if tscore > float(hypedsearch_score):
                    w.write("For " + filename + " DLQTLAL-NAAR scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore > float(hypedsearch_score):
                    w.write("For " + filename + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore > float(hypedsearch_score):
                    w.write("For " + filename + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALWSRM", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALWSRM", ppm_tolerance)
                if tscore > float(hypedsearch_score):
                    w.write("For " + filename + " DLQTLAL-WSRM scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
    
    w.write("\nWhen Hypedsearch found a target hybrid: \n \n")
            
    for file in spectra_files:
        filename = file.split("/")[-1]
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]
        prev_spec_id = -1

        hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
        for result in hypedsearch_answers:
            spec_id = result[0]
            hypedsearch_score = result[1]
            seq = result[2]
            if Multiset(seq) == Multiset("DLQTLAL-EVE") and spec_id != prev_spec_id:
                prev_spec_id = spec_id
                w.write("For " + filename + " Hypedsearch found DLQTLAL-EVE on spectrum: " + spec_id + " with score: " + hypedsearch_score + "\n")
            elif Multiset(seq) == Multiset("DLQTLAL-NAAR") and spec_id != prev_spec_id:
                prev_spec_id = spec_id
                w.write("For " + filename + " Hypedsearch found DLQTLAL-NAAR on spectrum: " + spec_id + " with score: " + hypedsearch_score + "\n")
            elif Multiset(seq) == Multiset("DLQTLAL-WSRM") and spec_id != prev_spec_id:
                prev_spec_id = spec_id
                w.write("For " + filename + " Hypedsearch found DLQTLAL-WSRM on spectrum: " + spec_id + " with score: " + hypedsearch_score + "\n")
            elif Multiset(seq) == Multiset("DPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                prev_spec_id = spec_id
                w.write("For " + filename+ " Hypedsearch found DPQVAQLELGG-EVEDPQVAQLELGGGPGAG on spectrum: " + spec_id + " with score: " + hypedsearch_score + "\n")
            elif Multiset(seq) == Multiset("EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                prev_spec_id = spec_id
                w.write("For " + filename + " Hypedsearch found EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG on spectrum: " + spec_id + " with score: " + hypedsearch_score + "\n")

    w.write("\nWhen Hypedsearch scored higher than Comet but target hybrid scored >8: \n \n")
    for file in spectra_files:
        filename = file.split("/")[-1]
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]
        prev_spec_id = -1

        hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
        for result in hypedsearch_answers:
            spec_id = result[0]
            hypedsearch_score = result[1]
            seq = result[2]
            hybrid_eval = result[3]
            input_spectrum = spectra[int(spec_id)]
            prec_tolerance = ppm_to_da(input_spectrum.precursor_mass, prec_tol)
            
            if abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALEVE", ppm_tolerance)
                if tscore < float(hypedsearch_score) and tscore > 8 and Multiset(seq) != Multiset("DLQTLAL-EVE") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-EVE scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALNAAR", ppm_tolerance)
                if tscore < float(hypedsearch_score) and tscore > 8 and Multiset(seq) != Multiset("DLQTLAL-NAAR") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-NAAR scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore < float(hypedsearch_score) and tscore > 8 and Multiset(seq) != Multiset("DPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore < float(hypedsearch_score) and tscore > 8 and Multiset(seq) != Multiset("EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALWSRM", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALWSRM", ppm_tolerance)
                if tscore < float(hypedsearch_score) and tscore > 8 and Multiset(seq) != Multiset("DLQTLAL-WSRM") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-WSRM scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")

              
    w.write("\nWhen Hypedsearch found something else with the same score as a target hybrid: \n \n")
            
    for file in spectra_files:
        filename = file.split("/")[-1]
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]
        prev_spec_id = -1

        hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
        for result in hypedsearch_answers:
            spec_id = result[0]
            hypedsearch_score = result[1]
            seq = result[2]
            hybrid_eval = result[3]
            input_spectrum = spectra[int(spec_id)]
            prec_tolerance = ppm_to_da(input_spectrum.precursor_mass, prec_tol)
            
            if abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALEVE", ppm_tolerance)
                if tscore == float(hypedsearch_score) and Multiset(seq) != Multiset("DLQTLA-LEVE") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-EVE scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALNAAR", ppm_tolerance)
                if tscore == float(hypedsearch_score) and Multiset(seq) != Multiset("DLQTLAL-NAAR") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-NAAR scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore == float(hypedsearch_score) and Multiset(seq) != Multiset("DPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                if tscore == float(hypedsearch_score) and Multiset(seq) != Multiset("EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALWSRM", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance: #and spec_id != prev_spec_id
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALWSRM", ppm_tolerance)
                if tscore == float(hypedsearch_score) and Multiset(seq) != Multiset("DLQTLAL-WSRM") and spec_id != prev_spec_id:
                    prev_spec_id = spec_id
                    w.write("For " + filename + " DLQTLAL-WSRM scored: " + str(tscore) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
                            
    w.write("\n Where all of the target hybrids could even appear \n \n")
    for file in spectra_files:
        filename = file.split("/")[-1]
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]
        
        hypedsearch_answers = get_hypedsearch_results(output_dict[int(file_id)])
        for input_spectrum in spectra:
            prec_tolerance = ppm_to_da(input_spectrum.precursor_mass, prec_tol)
            if abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALEVE", ppm_tolerance)
                w.write("For " + filename + " DLQTLAL-EVE is a match for spectrum: " + str(input_spectrum.num) + "." + "\tIt scores: " + str(tscore) + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALNAAR", ppm_tolerance)
                w.write("For " + filename + " DLQTLAL-NAAR is a match for spectrum: " + str(input_spectrum.num) + "." + "\tIt scores: " + str(tscore) + "\n")
            if abs(gen_spectra.get_precursor("DLQTLALWSRM", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DLQTLALWSRM", ppm_tolerance)
                w.write("For " + filename + " DLQTLAL-WSRM is a match for spectrum: " + str(input_spectrum.num) + "." + "\tIt scores: " + str(tscore) + "\n")
            if abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                w.write("For " + filename + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG is a match for spectrum: " + str(input_spectrum.num) + "." + "\tIt scores: " + str(tscore) + "\n")
            if abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tolerance:
                tscore = get_target_score(input_spectrum.mz_values, "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
                w.write("For " + filename + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG is a match for spectrum: " + str(input_spectrum.num) + "." + "\tIt scores: " + str(tscore) + "\n")
