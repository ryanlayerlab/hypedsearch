from preprocessing import preprocessing_utils, clustering, merge_search
from scoring import scoring
import gen_spectra
from main import get_spectra_files
from utils import ppm_to_da
import database
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS
from sqlite import database_file
import numpy as np
import matplotlib.pyplot as plt
import random

#These are user inputted parameters and their typical values
ppm_tolerance = 20 
peak_filter = 100 #given 25 peaks as information, some other software uses 50-100
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25 
make_new_db = False

def first_pass_truth_set(filepath):
    correct_sequences = []
    with open(filepath, 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split(';')
                correct_sequences.append(split_line[9])

    return correct_sequences

def find_overlaps(sequence, ppm_tol, input_masses):
    total_score = 0
    overlap_masses = list()
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
            overlap_masses.append(observed)
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return(total_score, overlap_masses)

#Set your filepaths to the database and the spectra folder
prot_path = '/home/karo9276/HypedSearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(prot_path)

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
specmill_seqs = first_pass_truth_set(truth_set_path)

#See which specmill seqs have len < 25
ctr = 0
for seq in specmill_seqs:
    if(len(seq) >= 25):
        ctr = ctr + 1

#In this dataset, these are known hybrids (added manually)
#hybrid_seqs = specmill_seqs[4:11]
#everything else in the dataset is a natural 

# #container stores result of if a spectrum is hybrid or not 
# is_hybrid = list()

# #container to store number of precursor matches for each spectrum
# b_precursor_matches = list()
# y_precursor_matches = list()

# #This loop checks each spectrum to determine if it is a hybrid peptide
# for i,spectrum in enumerate(spectra):
#     #This matches every mz in a spectrum with a list of kmers it can match to. Format is (m/z, location_start, location_end, ion, charge, parent_protein)
#     b_prec, y_prec = gen_spectra.convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
#     matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, b_prec, y_prec)

#     #FIRST ROUTE TO CHECK IF HYBRID - CHECK PRECURSOR WEIGHT MATCHES 
#     #Get precursor mass
#     b_precursor, y_precursor = list(matched_masses_b.keys())[-2], list(matched_masses_y.keys())[-1]

#     #Getting everything that matched the precursor weights
#     all_b_hits, all_y_hits = matched_masses_b[b_precursor], matched_masses_y[y_precursor]

#     #save number of precursor match hits
#     b_precursor_matches.append(len(all_b_hits))
#     y_precursor_matches.append(len(all_y_hits))
    
#     #Loop over all precursor weight hits and score each one 
#     b_scores, y_scores = list(), list()
#     seq_set = set()
#     for hit in all_b_hits: 
#         #Obtain sequence from hit 
#         seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
#         seq_set.add(seq)
#         #Score the hit - score = number of peaks matched (out of 25)
#         b_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
#     for hit in all_y_hits: 
#         seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
#         seq_set.add(seq)
#         y_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    
#     #If SpecMill sequence is found in our precursor match sequences anywhere, not hybrid? 
#     #What if seq_set is empty?
#     count = 0
#     if(specmill_seqs[i] in seq_set):
#         count = count + 1  

#     if(count == 0 ): 
#         is_hybrid.append(1)
#     else: 
#         is_hybrid.append(0)
        
#     #Perhaps if the mean scores are below some threshold, they are hybrid? 
#         #Need to truth this with some instances 

#     #SECOND ROUTE TO CHECK IF HYBRID - CHECK HIGHEST MASS MATCH 
#     #Best scoring, highest mass match (with charge accounted for...)
#     #How does it compare to the precursor weight? 
#     #If there is a mass close to the precursor with a 'high' score, it may not be hybrid 

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

known_hybrid_seqs = specmill_seqs[4:11]
known_natural_seqs = specmill_seqs[:4] + specmill_seqs[12:]

def rand_jitter(arr):
    stdev = .01 * (max(arr) - min(arr))
    random_list = list()
    for element in arr:
        random_num = random.randrange(0, 1000)
        random_list.append(random_num)
    return random_list 

#Hybrids 
all_overlap_masses_hyb = list()
all_yaxis = list()
scores_hyb = list() 
for i,hybrid in enumerate(known_hybrids):
    hybrid_seq = known_hybrid_seqs[i]
    score, overlap_masses = find_overlaps(hybrid_seq, ppm_tolerance, hybrid.mz_values)
    [all_overlap_masses_hyb.append(x) for x in overlap_masses]
    y_axis = rand_jitter(overlap_masses)
    [all_yaxis.append(x) for x in y_axis]
    [scores_hyb.append(score) for x in overlap_masses]
    
#Naturals 
all_overlap_masses_nat = list()
all_yaxis_nat = list()
scores_nat = list()
for i,natural in enumerate(known_naturals):
    natural_seq = known_natural_seqs[i]
    score, overlap_masses = find_overlaps(natural_seq, ppm_tolerance, natural.mz_values)
    [all_overlap_masses_nat.append(x) for x in overlap_masses]
    if len(overlap_masses) != 0: 
        y_axis = rand_jitter(overlap_masses)
        [all_yaxis_nat.append(x) for x in y_axis]
    [scores_nat.append(score) for x in overlap_masses]
    
plt1, ax1 = plt.subplots() #color can be red, blue, green, etc
ax1.scatter(all_overlap_masses_nat, all_yaxis_nat, color = 'b', alpha=.1, label='natural')
ax1.scatter(all_overlap_masses_hyb, all_yaxis, color = 'r', alpha=.5, label='hybrid')
plt.xlabel("Mz")
plt.ylabel("random jitter")
plt.title("Mass hits naturals vs hybrids")
plt.legend()
plt.savefig("mass_hits_norm")