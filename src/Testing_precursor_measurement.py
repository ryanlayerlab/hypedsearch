import database
from preprocessing import preprocessing_utils
from main import get_spectra_files
from utils import ppm_to_da
import matplotlib.pyplot as plt
from sqlite import database_file
from gen_spectra import convert_precursor_to_ion, get_precursor, gen_spectrum
import matplotlib.pyplot as plt

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 10

prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(prot_path)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/NOD2_E3'
spectra_files = get_spectra_files(spectra_path)
spectra, _ = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

def truth_set(filepath):
    correct_sequences = []
    with open(filepath, 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split(';')
                correct_sequences.append(split_line[9])
                
    return correct_sequences

specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/NOD2_E3_results.ssv'
correct_sequences = truth_set(specmill_path)

# def find_hit(protein_list, correct_seq):
#     for i, prot in enumerate(protein_list):
#         prot_seq = prot[1]
#         if correct_seq in prot_seq:
#             pid = i
#             break
#     for i in range(0,len(prot_seq)-len(correct_seq)):
#         if prot_seq[i:len(correct_seq)] == correct_seq:
#             return i
#     return -1 #Code should never go here

def get_scores(top_scores):
    scores = []
    for top in top_scores:
        scores.append(top[1])
        
    return scores
    
def test_existence(prec_hits, protein_list, specmill):
    for hit in prec_hits:
        start, end = hit[1], hit[2]
        pid = hit[5]
        prot_seq = protein_list[pid][1]
        hit_seq = prot_seq[start:end]
        if hit_seq == specmill:
            return True, hit
        
    return False, []

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

def test_scores(target_hit, all_hits, protein_list, obs_prec, prec_charge, input_masses, ppm_tol):
    scored = []
    for hit in all_hits:
        seq = find_sequence(hit, protein_list)
        distscore = distscoring(seq, obs_prec, prec_charge)
        aligned_score = overlap_scoring(seq, ppm_tol, input_masses)
        scored.append((aligned_score, distscore, seq, target_hit))
    
    scored = sorted(scored, key = lambda x: (x[0], x[1]), reverse=True)
    top_hit = scored[0]
    next_hit = scored[1]
    counter = 2
    while arePermutation(top_hit[2], next_hit[2]):
        next_hit = scored[counter]
        counter = counter + 1
    if scored[0][3] == target_hit:
        return True, top_hit[0] - next_hit[0], top_hit[0]
    else:
        return False, 0, 0

#See what we can identify from just the precursor mass
fails, finds = [], []
top_score, not_top, best_scores = [], [], []
for j,spectrum in enumerate(spectra):
    prec = spectrum.precursor_mass
    dbf = database_file(10, False)
    converted_b, converted_y = convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    conv_prec_tol = ppm_to_da(converted_b, prec_tol)
    prec_hits = dbf.query_mass(converted_b, conv_prec_tol)[0]
    existence, solution = test_existence(prec_hits, proteins.proteins, correct_sequences[j])
    if existence == False:
        fails.append(j)
        
    else:
        finds.append(j)
        score, distance, best = test_scores(solution, prec_hits, proteins.proteins, spectrum.precursor_mass, spectrum.precursor_charge, spectrum.mz_values, ppm_tolerance)
        if score == True:
            top_score.append((j, distance))
            best_scores.append((j, best))
        else:
            not_top.append(j)

print("There were", len(finds), "spectra where we found the specmill sequence from precursor")
print("Of these,", len(top_score), "spectra were the highest scoring among other precursor matches with an average distance of the top score to second best of,", sum(get_scores(top_score))/len(top_score))
print("There were", len(not_top), "spectra where the precursor was not the top scoring precursor choice")
print("There were", len(fails), "spectra where we didn't find the specmill sequence from precursor")
with open("fails_and_finds.txt", 'w') as f:
    f.write("Fails:\n")
    [f.write(str(x) + ', ') for x in fails]
    f.write("\n")
    f.write("Finds:\n")
    [f.write(str(x) + ', ') for x in finds]

item, score = [], []
for top in top_score:
    item.append(top[0])
    score.append(top[1])
sorted_score = sorted(score, reverse=True)
indices = []
[indices.append(x) for x in range(0,len(sorted_score))]
plt1, ax1 = plt.subplots()
ax1.scatter(indices, sorted_score, sizes = [10], alpha = .5)
plt.xlabel("Index")
plt.ylabel("Distance between top score and second best score")
plt.title("Distance between top score and second best score for each spectra")
plt.savefig("Distribution of top scores vs second highest")

spec, top = [], []
for score in best_scores:
    spec.append(score[0])
    top.append(score[1])
plt2, ax2 = plt.subplots()
ax2.scatter(spec, top, sizes = [10], alpha=.5)
plt.xlabel("Spectrum Index")
plt.ylabel("Score")
plt.title("The value of the top score for each spectrum we found from just the precursor")
plt.savefig("Spectrum Idx")