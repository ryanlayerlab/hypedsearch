import database
import os
import collections
import numpy as np
from preprocessing import preprocessing_utils, clustering, merge_search
import gen_spectra
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import modified_match_masses
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sqlite import database_file
from gen_spectra import get_raw_mass
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25
digest = True

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/Lab_Data'
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
        prev_id=-1
        for line in f:
            #spec_id	Hybrid?	seq	total_score	precursor_score	left_parent	right_parent	left_score	right_score	25_extended_seq
            A = line.split("\t")
            spec_id = A[0]
            if spec_id != prev_id:
                prev_id = spec_id
                score = A[3]
                seq = A[2]
                eval = A[1]
                result.append((spec_id, score, seq, eval))
    return result

def get_target_scores(input_masses, target1, target2, target3, target4, ppm_tol):
    target1_score = overlap_scoring(target1, input_masses, ppm_tol, )
    target2_score = overlap_scoring(target2, input_masses, ppm_tol)
    target3_score = overlap_scoring(target3, input_masses, ppm_tol)
    target4_score = overlap_scoring(target4, input_masses, ppm_tol)
    return target1_score, target2_score, target3_score, target4_score

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

output_dict = get_output_files("/home/naco3124/jaime_hypedsearch/hypedsearch/data/output")

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
    for file in spectra_files:
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
            prec_tol = ppm_to_da(input_spectrum.precursor_mass, ppm_tolerance)
            #score all the good ones and see if it beats hypedsearch score
            score1, score2, score3, score4 = get_target_scores(input_spectrum.mz_values, "DLQTLALEVE", "DLQTLALNAAR",  "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
            if score1 > int(hypedsearch_score) and abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DLQTLALEVE scored: " + str(score1) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score2 > int(hypedsearch_score) and abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DLQTLALNAAR scored: " + str(score2) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score3 > int(hypedsearch_score) and abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(score3) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score4 > int(hypedsearch_score) and abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(score3) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
                    
    w.write("\nNow looking at the case of equality \n \n")
            
    for file in spectra_files:
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
            prec_tol = ppm_to_da(input_spectrum.precursor_mass, ppm_tolerance)
            #score all the good ones and see if it beats hypedsearch score
            score1, score2, score3, score4 = get_target_scores(input_spectrum.mz_values, "DLQTLALEVE", "DLQTLALNAAR",  "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
            if score1 == int(hypedsearch_score) and abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DLQTLALEVE scored: " + str(score1) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score2 == int(hypedsearch_score) and abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DLQTLALNAAR scored: " + str(score2) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score3 == int(hypedsearch_score) and abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(score3) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")
            if score4 == int(hypedsearch_score) and abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + file_id + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored: " + str(score4) + " on spectrum: " + spec_id + " while Hypedsearch returned: " + seq + " as a: " + hybrid_eval + " with score: " + hypedsearch_score + "\n")

    w.write("\nNow looking at the case of things Spectrummill scored >8 \n \n")
    for file in spectra_files:
        spectra = preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        filesplit = file.split("Fxn")[1]
        file_id = filesplit.split(".")[0]
        for spec_id, input_spectrum in enumerate(spectra):
            score1, score2, score3, score4 = get_target_scores(input_spectrum.mz_values, "DLQTLALEVE", "DLQTLALNAAR",  "DPQVAQLELGGEVEDPQVAQLELGGGPGAG", "EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", ppm_tolerance)
            if score1 > 8 and abs(gen_spectra.get_precursor("DLQTLALEVE", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + str(file_id) + " DLQTLAL-EVE scored " + str(score1) + " for spectrum: " + str(spec_id) + "\n")
            if score2 > 8 and abs(gen_spectra.get_precursor("DLQTLALNAAR", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + str(file_id) + " DLQTLAL-NAAR scored " + str(score2) + " for spectrum: " + str(spec_id) + "\n")
            if score3 > 8 and abs(gen_spectra.get_precursor("DPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + str(file_id) + " DPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored " + str(score3) + " for spectrum: " + str(spec_id) + "\n")
            if score4 > 8 and abs(gen_spectra.get_precursor("EVEDPQVAQLELGGEVEDPQVAQLELGGGPGAG", input_spectrum.precursor_charge) - input_spectrum.precursor_mass) < prec_tol:
                w.write("For BMEM: " + str(file_id) + " EVEDPQVAQLELGG-EVEDPQVAQLELGGGPGAG scored " + str(score4) + " for spectrum: " + str(spec_id) + "\n")