from preprocessing import preprocessing_utils, clustering, merge_search
from scoring import scoring
import gen_spectra
from main import get_spectra_files
from lookups.utils import ppm_to_da
import computational_pipeline.database_generator as database_generator
from lookups.constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS
from sqlite import database_file
import numpy as np
import matplotlib.pyplot as plt
import random

#These are user inputted parameters and their typical values
ppm_tolerance = 20 
peak_filter = 50 #given 25 peaks as information, some other software uses 50-100
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25 
make_new_db = False

def first_pass_truth_set(filepath):
    specmill_nat, specmill_hyb = [], []
    with open(filepath, 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split(';')
                if 'HYBRID' in split_line[15]:
                    info = split_line[15].split(' ')
                    info_split = info[3].split('-')
                    for i in range(0,len(split_line[9])):
                        if split_line[9][0:i] in info_split[0]:
                            hyb_seq = (split_line[9][0:i])
                        else:
                            hyb_junc = hyb_seq + '-' + split_line[9][(i-1):]
                            specmill_hyb.append(hyb_junc)
                            break
                else:
                    specmill_nat.append(split_line[9])

    return specmill_nat, specmill_hyb

def find_overlaps(masses, ppm_tol, input_masses):
    total_score = 0
    overlap_masses = list()
    masses = sorted(masses['spectrum'])
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
            overlap_masses.append(observed)
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return((total_score/len(input_masses)), overlap_masses)

#Set your filepaths to the database and the spectra folder
prot_path = '/home/karo9276/HypedSearch/hypedsearch/data/database/sample_database.fasta'
proteins = database_generator.build_database(prot_path)

dbf = database_file(max_pep_len, make_new_db)

if make_new_db:
    kv_prots = [(k, v) for k, v in proteins.proteins]    
    merge_search.modified_make_database_set(kv_prots, max_pep_len, dbf)

spectra_path = '/home/karo9276/HypedSearch/hypedsearch/data/spectra/NOD2_E3'
spectra_file = get_spectra_files(spectra_path)[0]

#Loads in the spectra as a list of spectrum objects
spectra = preprocessing_utils.load_spectra(spectra_file, ppm_tolerance, peak_filter, relative_abundance_filter)

#Loading in the truth set from SpectraMill
truth_set_path = '/home/karo9276/HypedSearch/hypedsearch/data/truth_table/NOD2_E3_results.ssv'
specmill_natural_seqs, specmill_hybrid_seqs = first_pass_truth_set(truth_set_path)

#In this dataset, these are known hybrids (added manually)
#everything else in the dataset is a natural 

# #container stores result of if a spectrum is hybrid or not 
is_hybrid = list()
is_hybrid_seq = list()

# #container to store number of precursor matches for each spectrum
b_precursor_matches = list()
y_precursor_matches = list()

""" # #This loop checks each spectrum to determine if it is a hybrid peptide
for i,spectrum in enumerate(spectra):
    #This matches every mz in a spectrum with a list of kmers it can match to. Format is (m/z, location_start, location_end, ion, charge, parent_protein)
    b_prec, y_prec = gen_spectra.convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, b_prec, y_prec)

    #FIRST ROUTE TO CHECK IF HYBRID - CHECK PRECURSOR WEIGHT MATCHES 
#    #Get precursor mass
    b_precursor, y_precursor = list(matched_masses_b.keys())[-2], list(matched_masses_y.keys())[-1]

#   #Getting everything that matched the precursor weights
    all_b_hits, all_y_hits = matched_masses_b[b_precursor], matched_masses_y[y_precursor]

#    #save number of precursor match hits
    b_precursor_matches.append(len(all_b_hits))
    y_precursor_matches.append(len(all_y_hits))
    
#    #Loop over all precursor weight hits and score each one 
    b_scores, y_scores = list(), list()
    seq_set = set()
    for hit in all_b_hits: 
#        #Obtain sequence from hit 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        seq_set.add(seq)
        #Score the hit - score = number of peaks matched (out of 25)
        #b_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    for hit in all_y_hits: 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        seq_set.add(seq)
        #y_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
         #If SpecMill sequence is found in our precursor match sequences anywhere, not hybrid? 
    #What if seq_set is empty?
    count = 0
    if(specmill_seqs[i] in seq_set):
        count = count + 1  

    if(count == 0 ): 
        is_hybrid.append(1)
    else: 
        is_hybrid.append(0)
        
    if(len(seq_set) == 0 & len(specmill_seqs[i]) >= 25):
        is_hybrid_seq.append(1)
    else: 
        is_hybrid_seq.append(0)
#     #Perhaps if the mean scores are below some threshold, they are hybrid? 
#         #Need to truth this with some instances 

#     #SECOND ROUTE TO CHECK IF HYBRID - CHECK HIGHEST MASS MATCH 
#     #Best scoring, highest mass match (with charge accounted for...)
#     #How does it compare to the precursor weight? 
#     #If there is a mass close to the precursor with a 'high' score, it may not be hybrid 
 """
# print("The number of potential hybrids is: ")
# print(sum(is_hybrid))
# print(" out of: ")
# print(len(is_hybrid))

# with open("hybrid_indices.txt", "w") as h:
#     for i,entry in enumerate(is_hybrid):
#         if entry != 0:
#             h.write(str(i) + "\n")

# missed_seqs = list()
# with open("hybrid_indices.txt", "r") as f:
#     for line in f:
#         if(len(specmill_seqs[int(line)]) <= 25):
#             missed_seqs.append(line)
            
# with open("missed_indices.txt", "w") as h:
#     for entry in missed_seqs:
#         h.write(entry)

#DETERMINE HOW BIG ARE GOOD HITS IN HYBRIDS VS. NATURALS 
known_hybrids = spectra[4:11]
known_naturals = spectra[:4] + spectra[12:]

known_hybrid_seqs = specmill_hybrid_seqs
known_natural_seqs = specmill_natural_seqs

#Hybrids 
all_overlap_masses_hyb = list()
biggest_b_hit_hyb = list()
biggest_y_hit_hyb = list()
all_yaxis = list()
scores_hyb = list() 
for i,hybrid in enumerate(known_hybrids):
    hybrid_seq = known_hybrid_seqs[i]
    #split to left and right pieces 
    hybrid_seq_split = hybrid_seq.split('-')
    #if len(hybrid_seq) <= 25:
    spectrum_b_1 = gen_spectra.generate_spectrum(hybrid_seq_split[0], 1, 'b')
    spectrum_b_2 = gen_spectra.generate_spectrum(hybrid_seq_split[0], 2, 'b')
    spectrum_y_1 = gen_spectra.generate_spectrum(hybrid_seq_split[1], 1, 'y')
    spectrum_y_2 = gen_spectra.generate_spectrum(hybrid_seq_split[1], 2, 'y')
    score_b1, overlap_masses_b1 = find_overlaps(spectrum_b_1, ppm_tolerance, hybrid.mz_values)
    score_b2, overlap_masses_b2 = find_overlaps(spectrum_b_2, ppm_tolerance, hybrid.mz_values)
    score_y1, overlap_masses_y1 = find_overlaps(spectrum_y_1, ppm_tolerance, hybrid.mz_values)
    score_y2, overlap_masses_y2 = find_overlaps(spectrum_y_2, ppm_tolerance, hybrid.mz_values)
        
    converted_masses_b1 = list()
    converted_masses_b2 = list()
    converted_masses_y1 = list()
    converted_masses_y2 = list()
    for mass in overlap_masses_b1:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 0, 1, hybrid.precursor_charge)
        converted_masses_b1.append(converted_mass)
    for mass in overlap_masses_b2:
            converted_mass = gen_spectra.convert_ion_to_precursor(mass, 0, 2, hybrid.precursor_charge)
            converted_masses_b2.append(converted_mass)
    for mass in overlap_masses_y1:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 1, 1, hybrid.precursor_charge)
        converted_masses_y1.append(converted_mass)
    for mass in overlap_masses_y2:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 1, 2, hybrid.precursor_charge)
        converted_masses_y2.append(converted_mass)
        
    #Get the biggest hit for b and y ions 
    if len(converted_masses_b1) == 0 and len(converted_masses_b2) == 0:
        biggest_b_hit = 0
    else:
        biggest_b_hit = max(converted_masses_b1 + converted_masses_b2) #biggest b hit 
    if len(converted_masses_y1) == 0 and len(converted_masses_y2) == 0:
        biggest_y_hit = 0
    else: 
        biggest_y_hit = max(converted_masses_y1 + converted_masses_y2)   
        
    #Divide by the precursor mass 
    hyb_prec = hybrid.precursor_mass
    hyb_charge = hybrid.precursor_charge
    b_hit_prop = biggest_b_hit/hyb_prec
    y_hit_prop = biggest_y_hit/hyb_prec
    biggest_b_hit_hyb.append(b_hit_prop)
    biggest_y_hit_hyb.append(y_hit_prop)
    scores_hyb.append(score_b1 + score_b2 + score_y1 + score_y2)
             
#Naturals 
all_overlap_masses_nat = list()
biggest_b_hit_nat = list()
biggest_y_hit_nat = list()
all_yaxis_nat = list()
scores_nat = list()
for i,natural in enumerate(known_naturals):
    natural_seq = known_natural_seqs[i]
    #if len(natural_seq) <= 25:
    spectrum_b_1 = gen_spectra.generate_spectrum(natural_seq, 1, 'b')
    spectrum_b_2 = gen_spectra.generate_spectrum(natural_seq, 2, 'b')
    spectrum_y_1 = gen_spectra.generate_spectrum(natural_seq, 1, 'y')
    spectrum_y_2 = gen_spectra.generate_spectrum(natural_seq, 2, 'y')
    score_b1, overlap_masses_b1 = find_overlaps(spectrum_b_1, ppm_tolerance, natural.mz_values)
    score_b2, overlap_masses_b2 = find_overlaps(spectrum_b_2, ppm_tolerance, natural.mz_values)
    score_y1, overlap_masses_y1 = find_overlaps(spectrum_y_1, ppm_tolerance, natural.mz_values)
    score_y2, overlap_masses_y2 = find_overlaps(spectrum_y_2, ppm_tolerance, natural.mz_values)
        
    converted_masses_b1 = list()
    converted_masses_b2 = list()
    converted_masses_y1 = list()
    converted_masses_y2 = list()
    for mass in overlap_masses_b1:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 0, 1, natural.precursor_charge)
        converted_masses_b1.append(converted_mass)
    for mass in overlap_masses_b2:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 0, 2, natural.precursor_charge)
        converted_masses_b2.append(converted_mass)
    for mass in overlap_masses_y1:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 1, 1, natural.precursor_charge)
        converted_masses_y1.append(converted_mass)
    for mass in overlap_masses_y2:
        converted_mass = gen_spectra.convert_ion_to_precursor(mass, 1, 2, natural.precursor_charge)
        converted_masses_y2.append(converted_mass)
        
    #Get the biggest hit for b and y ions 
    if len(converted_masses_b1) == 0 and len(converted_masses_b2) == 0:
        biggest_b_hit = 0
    else:
        biggest_b_hit = max(converted_masses_b1 + converted_masses_b2) #biggest b hit 
    if len(converted_masses_y1) == 0 and len(converted_masses_y2) == 0:
        biggest_y_hit = 0
    else: 
        biggest_y_hit = max(converted_masses_y1 + converted_masses_y2)   
     
    #Divide by the precursor mass 
    nat_prec = natural.precursor_mass
    nat_charge = natural.precursor_charge
    biggest_b_hit_nat.append(biggest_b_hit/nat_prec)
    biggest_y_hit_nat.append(biggest_y_hit/nat_prec)
    scores_nat.append(score_b1 + score_b2 + score_y1 + score_y2)

#make axes have the same scale, make plots a bit larger, add mean lines 
plt1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8,6))
ax1.scatter(biggest_b_hit_nat, scores_nat, color='g', label='natural, b', alpha=0.6)
ax1.axvline(sum(biggest_b_hit_nat)/len(biggest_b_hit_nat))
ax1.set_xlim(0,1)
ax1.set_ylim(0,0.6)
ax3.scatter(biggest_y_hit_nat, scores_nat, color='g', label='natural, y', alpha=0.6)
ax3.axvline(sum(biggest_y_hit_nat)/len(biggest_y_hit_nat))
ax3.set_xlim(0,1)
ax3.set_ylim(0,0.6)
ax2.scatter(biggest_b_hit_hyb, scores_hyb, color='r', label='hybrid, b')
ax2.axvline(sum(biggest_b_hit_hyb)/len(biggest_b_hit_hyb))
ax2.set_xlim(0,1)
ax2.set_ylim(0,0.6)
ax4.scatter(biggest_y_hit_hyb, scores_hyb, color='r', label='hybrid, y')
ax4.axvline(sum(biggest_y_hit_hyb)/len(biggest_y_hit_hyb))
ax4.set_xlim(0,1)
ax4.set_ylim(0,0.6)
plt.xlabel("Max. mapped mass / precursor mass")
plt.ylabel("Score")
plt.suptitle("Maximum mapped mass as a proportion of the total peptide mass")
plt.savefig("mass_hit_score_normalized")

plt2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(8,6))
ax1.hist(biggest_b_hit_nat, bins=30, color='g', label='natural, b', alpha=0.6)
ax1.axvline(sum(biggest_b_hit_nat)/len(biggest_b_hit_nat))
ax1.set_xlim(0,1)
ax3.hist(biggest_y_hit_nat, bins=30, color='g', label='natural, y', alpha=0.6)
ax3.axvline(sum(biggest_y_hit_nat)/len(biggest_y_hit_nat))
ax3.set_xlim(0,1)
ax2.hist(biggest_b_hit_hyb, bins=30, color='r', label='hybrid, b', alpha=0.6)
ax2.axvline(sum(biggest_b_hit_hyb)/len(biggest_b_hit_hyb))
ax2.set_xlim(0,1)
ax4.hist(biggest_y_hit_hyb, bins=30, color='r', label='hybrid, y', alpha=0.6)
ax4.axvline(sum(biggest_y_hit_hyb)/len(biggest_y_hit_hyb))
ax4.set_xlim(0,1)
plt.xlabel("Max. mapped mass / precursor mass")
plt.ylabel("Freq.")
plt.suptitle("Maximum mapped mass as a proportion of the total peptide mass")
plt.savefig("mass_hit_score_normalized_hist")

plt3, (ax1,ax2) = plt.subplots(2, figsize=(5,4))
ax1.hist(scores_nat, bins=30, color='g', label='natural', alpha=0.6)
ax1.set_xlim(0,0.5)
ax2.hist(scores_hyb, bins=30, color='r', label='hybrid', alpha=0.6)
ax2.set_xlim(0,0.5)
plt.xlabel("Proportion of theoretical peaks matched")
plt.ylabel("Freq.")
plt.suptitle("Proportion of peaks matched")
plt.legend()
plt.savefig("score_prop")
       
#plt1, (ax1, ax2) = plt.subplots(2)
#ax1.scatter(range(0,len(precursor_score_hyb)), precursor_score_hyb, color='r', label='hybrid')
#ax2.scatter(range(0,len(precursor_score_nat)), precursor_score_nat, color='g', label='natural')
#plt.xlabel("Score")
#plt.ylabel("Frequency")
#plt.legend()
#plt.suptitle("Precursor scores for hybrid and naturals")
#plt.savefig("nat_hyb_scores")
    
#plt1, ax1 = plt.subplots() #color can be red, blue, green, etc
#ax1.scatter(all_overlap_masses_nat, all_yaxis_nat, color = 'b', alpha=.1, label='natural')
#ax1.scatter(all_overlap_masses_hyb, all_yaxis, color = 'r', alpha=.5, label='hybrid')
#plt.xlabel("Mz (/len(seq))")
#plt.ylabel("random jitter")
#plt.title("Mass hits naturals vs hybrids")
#plt.legend()
#plt.savefig("mass_hits_norm")

#plt2, ax1 = plt.subplots() #color can be red, blue, green, etc
#ax1.scatter(all_overlap_masses_nat, scores_nat, color = 'g', alpha=.1, label='natural')
#ax1.scatter(all_overlap_masses_hyb, scores_hyb, color = 'r', alpha=0.75, label='hybrid')
#plt.xlabel("mz")
#plt.ylabel("Proportion of peaks matched for each spectrum")
#plt.title("Proportion of peaks matched")
#plt.legend()
#plt.savefig("matched_peaks_size")

#plt3, ax1 = plt.subplots()
#ax1.hist(scores_nat, color='b', alpha=0.5, label='natural')
#ax1.hist(scores_hyb, color='r', alpha=0.5, label='hybrid')
#plt.xlabel("Score")
#plt.ylabel("Freq.")
#plt.legend()
#plt.savefig("score_hist")

#plt3, ax1 = plt.subplots()
#ax1.hist(all_overlap_masses_nat, color='b', alpha=0.5, label='natural')
#ax1.hist(all_overlap_masses_hyb, color='r', alpha=0.5, label='hybrid')
#plt.xlabel("Size")
#plt.ylabel("Freq.")
#plt.legend()
#plt.axvline(sum(all_overlap_masses_nat)/len(all_overlap_masses_nat))
#plt.savefig("size_hist")


#2-21 experiments
#percent_matched = list()
#garbage_matches = list()
#missing_matched = list()
#for i,spectrum in enumerate(spectra):
#    #Checking proportion of peaks matched for spectrummill seqs
#    specmill_seq = specmill_seqs[i]
#    score, _ = find_overlaps(specmill_seq, ppm_tolerance, spectrum.mz_values)
#    percent_matched.append(score/peak_filter)
#    missing_matched.append(score/(4*len(specmill_seq)))
#    
##    #Finding # garbage peaks 
#    b_prec, y_prec = gen_spectra.convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
#    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, b_prec, y_prec)
    #check how many things match to nothing 
#    count = 0
#    for mass in spectrum.mz_values:
#        if matched_masses_b[mass] == [] and matched_masses_y[mass] == []:
#            count = count + 1  
#    garbage_matches.append(count/len(spectrum.mz_values))
#    
#plotting the percent of peaks matched 
#plt4, ax1 = plt.subplots() #color can be red, blue, green, etc
#ax1.hist(percent_matched, color = 'g', alpha=0.5) #range(0,len(spectra))
#plt.xlabel("Proportion of peaks matched")
#plt.ylabel("Frequency")
#plt.title("Proportion of expected peaks matched for each spectrum")
#plt.savefig("peak_match")

#plotting proportion of missing expected sequences 
#plt4, ax1 = plt.subplots() #color can be red, blue, green, etc
#ax1.hist(missing_matched, color = 'g', alpha=0.5)
#plt.xlabel("Proportion of possible peaks found")
#plt.ylabel("Frequency")
#plt.title("Proportion of possible peaks found for each spectrum")
#plt.savefig("missing_match")

#Plotting the number of garbage peaks (didn't match to anything)
#plt5, ax1 = plt.subplots() #color can be red, blue, green, etc
#ax1.hist(garbage_matches, color = 'g', alpha=0.5)
#plt.xlabel("Proportion of 'noisy' peaks")
#plt.ylabel("Frequency")
#plt.title("Proportion of noisy peaks for each spectrum")
#plt.savefig("garb_peaks")