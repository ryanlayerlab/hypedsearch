from preprocessing import preprocessing_utils, clustering, merge_search
from scoring import scoring
import gen_spectra
from main import get_spectra_files
from utils import ppm_to_da
import database
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS
from sqlite import database_file
#import matplotlib.pyplot as plt

#These are user inputted parameters and their typical values
ppm_tolerance = 20 
peak_filter = 25 #given 25 peaks as information, some other software uses 50-100
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

#Set your filepaths to the database and the spectra folder
prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(prot_path)

dbf = database_file(max_pep_len, make_new_db)

if make_new_db:
    kv_prots = [(k, v) for k, v in proteins.proteins]    
    merge_search.modified_make_database_set(kv_prots, max_pep_len, dbf)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/NOD2_E3'
spectra_file = get_spectra_files(spectra_path)[0]

#Loads in the spectra as a list of spectrum objects
spectra = preprocessing_utils.load_spectra(spectra_file, ppm_tolerance, peak_filter, relative_abundance_filter)

#Loading in the truth set from SpectraMill
truth_set_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/NOD2_E3_results.ssv'
specmill_seqs = first_pass_truth_set(truth_set_path)

#See which specmill seqs have len < 25
ctr = 0
for seq in specmill_seqs:
    if(len(seq) >= 25):
        ctr = ctr + 1

#In this dataset, these are known hybrids (added manually)
#hybrid_seqs = specmill_seqs[4:11]
#everything else in the dataset is a natural 

#container stores result of if a spectrum is hybrid or not 
is_hybrid = list()

#This loop checks each spectrum to determine if it is a hybrid peptide
for i,spectrum in enumerate(spectra):
    #This matches every mz in a spectrum with a list of kmers it can match to. Format is (m/z, location_start, location_end, ion, charge, parent_protein)
    b_prec, y_prec = gen_spectra.convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, b_prec, y_prec)

    #FIRST ROUTE TO CHECK IF HYBRID - CHECK PRECURSOR WEIGHT MATCHES 
    #Get precursor mass
    b_precursor, y_precursor = list(matched_masses_b.keys())[-2], list(matched_masses_y.keys())[-1]

    #Getting everything that matched the precursor weights
    all_b_hits, all_y_hits = matched_masses_b[b_precursor], matched_masses_y[y_precursor]

    #Loop over all precursor weight hits and score each one 
    b_scores, y_scores = list(), list()
    seq_set = set()
    for hit in all_b_hits: 
        #Obtain sequence from hit 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        seq_set.add(seq)
        #Score the hit - score = number of peaks matched (out of 25)
        b_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    for hit in all_y_hits: 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        seq_set.add(seq)
        y_scores.append(overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    
    #If SpecMill sequence is found in our precursor match sequences anywhere, not hybrid? 
    #What if seq_set is empty?
    count = 0
    if(specmill_seqs[i] in seq_set):
        count = count + 1  

    if(count == 0 ): 
        is_hybrid.append(1)
    else: 
        is_hybrid.append(0)
        
    #Perhaps if the mean scores are below some threshold, they are hybrid? 
        #Need to truth this with some instances 

    #SECOND ROUTE TO CHECK IF HYBRID - CHECK HIGHEST MASS MATCH 
    #Best scoring, highest mass match (with charge accounted for...)
    #How does it compare to the precursor weight? 
    #If there is a mass close to the precursor with a 'high' score, it may not be hybrid 

print("The number of potential hybrids is: ")
print(sum(is_hybrid))
print(" out of: ")
print(len(is_hybrid))

with open("hybrid_indices.txt", "w") as h:
    for index in is_hybrid:
        if index != 0:
            h.write(str(index) + "\n")