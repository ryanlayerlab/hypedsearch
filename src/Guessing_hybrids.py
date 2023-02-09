from preprocessing import preprocessing_utils, clustering, merge_search
from scoring import scoring
from main import get_spectra_files
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

#Set your filepaths to the database and the spectra folder
prot_path = '/home/karo9276/Hybrid-Testing/Hybrid-Testing/data/database/sample_database.fasta'
proteins = database.build(prot_path)

dbf = database_file(max_pep_len, make_new_db)

if make_new_db:
    kv_prots = [(k, v) for k, v in proteins.proteins]    
    merge_search.modified_make_database_set(kv_prots, max_pep_len, dbf)

spectra_path = '/home/karo9276/Hybrid-Testing/Hybrid-Testing/data/spectra/NOD2_E3'
spectra_files = get_spectra_files(spectra_path)

#Loads in the spectra as a list of spectrum objects
spectra = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

#Loading in the truth set from SpectraMill
truth_set_path = '/home/karo9276/Hybrid-Testing/Hybrid-Testing/data/NOD2_E3_results.ssv'
specmill_seqs = preprocessing_utils.first_pass_truth_set(truth_set_path)

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
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, dbf)

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
        b_scores.append(scoring.overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    for hit in all_y_hits: 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        seq_set.add(seq)
        y_scores.append(scoring.overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    
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