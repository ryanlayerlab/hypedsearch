import database
import os
import collections
import numpy as np
from preprocessing import preprocessing_utils, clustering, merge_search
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import get_modified_match_masses
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

# prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
# proteins = database.build(prot_path)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/Lab_Data'
spectra_files = get_spectra_files(spectra_path)
spectra = []
for file in spectra_files:
    [spectra.append(x) for x in preprocessing_utils.load_spectra(file, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)]

def truth_set(filepath):
    correct_sequences = []
    with open(filepath, 'r') as truth_set:
        for line in truth_set:
            split_line = line.split('\t')
            correct_sequences.append(split_line[21])
            
    return correct_sequences

def get_score_map(hybrid_set, input_dict):
    score_dict = dict()
    for hybrid in hybrid_set:
        scores = input_dict[hybrid]
        for score in scores:
            if score not in score_dict.keys():
                score_dict[score] = []
            score_dict[score].append(hybrid)
    scores = [float(x) for x in score_dict.keys()]
    numbers = [len(score_dict[x]) for x in score_dict.keys()]
    return scores, numbers

def get_natives_and_hybrids(filepath):
    natives, hybrids = dict(), dict()
    with open(filepath, 'r') as truth_set:
        for i, line in enumerate(truth_set):
            split_line = line.split('\t')
            if split_line[21] == "DFLHAR":
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

if digest == False:
    specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/Adjusted_BMEM_searches.txt'
else:
    specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/Adjusted_Digest_BMEM_searches.txt'
correct_sequences = truth_set(specmill_path)
first_pass_specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/BMEM_searches.txt'
first_pass_specmill_seqs = truth_set(first_pass_specmill_path)
specmill_natives, native_scores, specmill_hybrids, hybrid_scores = get_natives_and_hybrids(specmill_path)
for seq in correct_sequences:
    if seq not in specmill_natives:
        if seq not in specmill_hybrids:
            print(seq)

# How many times did spectrumMill find the same thing?
# def hypedsearch_top_answers(output_files):
#     hhybrids, hnaturals = [], []
#     for file in output_files:
#         with open(file, 'r') as f:
#             prev_id = -1
#             for line in f:
#                 A = line.split('\t')
#                 if A[0] != prev_id:
#                     prev_score = A[3]
#                     prev_id = A[0]
#                 if A[3] == prev_score:
#                     prev_score = A[3]
#                     if A[1] == "Hybrid":
#                         hypedsearch_sequence = A[2].replace("-", "")
#                         hhybrids.append(hypedsearch_sequence)
#                     else:
#                         hypedsearch_sequence = A[2]
#                         hnaturals.append(hypedsearch_sequence)
#     return hnaturals, hhybrids

def same_permutation(a, b):
    d = collections.defaultdict(int)
    for x in a:
        d[x] += 1
    for x in b:
        d[x] -= 1
    return not any(d.itervalues())

def get_nonisomorphic_set(sequences):
    count = dict()
    new_seqs = [sequence[0]]
    for i, sequence in enumerate(sequences):
        for oseq in sequences[i:]:
            if same_permutation(sequence, oseq):
                sequences
                
def num_of_id_seqs(output_files):
    ncount, hcount = 0, 0
    for file in output_files:
        with open(file, 'r') as r:
            for line in r:
                A = line.split("\t")
                eval = A[1]
                if eval == "Hybrid":
                    hcount += 1
                else:
                    ncount += 1
    return ncount, hcount

def hypedsearch_all_answers(output_files):
    natives, hybrids = dict(), dict()
    for file in output_files:
        with open(file, 'r') as f:
            for line in f:
                A = line.split('\t')
                hypedsearch_sequence = A[2].replace("-", "")
                if A[1] == "Natural":                    
                    if hypedsearch_sequence not in natives.keys():
                        natives[hypedsearch_sequence] = []
                    natives[hypedsearch_sequence].append(A[3])
                else:                
                    if hypedsearch_sequence not in hybrids.keys():
                        hybrids[hypedsearch_sequence] = []
                    hybrids[hypedsearch_sequence].append(A[3])
                    
    return natives.keys(), natives, hybrids.keys(), hybrids

def test_isomorphic_seqs(seq1, seq2):
    #seqs are isomorphic if they are permutations or are I&L swapped
    seq1_set, seq2_set = set(), set()
    for char in seq1:
        seq1_set.add(char)
    for char in seq2:
        seq2_set.add(char)
    if len(seq1_set) == len(seq2_set):
        testperm = True
    else:
        testperm = False
        
    #Now to test if we can swap I & L
    lswap1 = seq1.replace("i", "l")
    lswap2 = seq2.replace("i", "l")
    if lswap1 == lswap2:
        testiso = True
    else:
        testiso = False
        
    if testperm and testiso:
        return True
    else:
        return False
    
def get_output_files(open_path):
    output_files = []
    for (root, _, filenames) in os.walk(open_path):
        for fname in filenames:
            output_files.append(os.path.join(root, fname))
    return output_files
    
open_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/output"
output_files = get_output_files(open_path)
hnatives, hnative_scores, hhybrids, hhybrid_scores = hypedsearch_all_answers(output_files)
hnat_num, hhyb_num = num_of_id_seqs(output_files)

def compare_with_SpecMill(hypedsearch_answer_set, correct_sequences, target_dict, hybrid_seqs):
    hcount = 0
    hcount_array = set()
    for answer in hypedsearch_answer_set:
        if answer in correct_sequences:
            if answer in hybrid_seqs:
                print("error")
            else:
                hcount = hcount + len(target_dict[answer])
                hcount_array.add(answer)
    return hcount, hcount_array

def compare_with_SpecMill_hybrids(hypedsearch_answer_set, correct_sequences, target_dict, native_seqs):
    hcount = 0
    hcount_array = set()
    for answer in hypedsearch_answer_set:
        if answer in correct_sequences:
            if answer in native_seqs:
                print("error")
            else:
                hcount = hcount + len(target_dict[answer])
                hcount_array.add(answer)
    return hcount, hcount_array

# def compare_unique_with_SpecMill(hypedsearch_answer_set, correct_sequences):
#     hcount = 0
#     for answer in hypedsearch_answer_set:
#         for specanswer in correct_sequences:
#             if test_isomorphic_seqs(answer, specanswer):
#                 hcount = hcount + 1
#     return hcount
    
hcountn, specmill_matched_natives = compare_with_SpecMill(hnatives, correct_sequences, native_scores, specmill_hybrids)
hcounth, specmill_new_hybrids = compare_with_SpecMill_hybrids(hhybrids, correct_sequences, hybrid_scores, specmill_natives)
# hucounth, _ = compare_with_SpecMill(huhybrids, correct_sequences)
# hucountn, _ = compare_with_SpecMill(hunaturals, correct_sequences)

def check_for_changed_mind(first_pass_truth_set, second_pass_truth_set):
    count = 0
    for seq in first_pass_truth_set:
        if seq not in second_pass_truth_set:
            count = count + 1
    return count

changednum = check_for_changed_mind(first_pass_specmill_seqs, correct_sequences)

scount = 0
nacount = 0
hycount = 0
for answer in correct_sequences:
    if answer in hnatives:
        nacount = nacount + 1
    elif answer in hhybrids:
        hycount = hycount + 1
scount = nacount + hycount
        
hcount = 0
for hybrid in hhybrids:
    hcount += len(hhybrid_scores[hybrid])
ncount = 0
for native in hnatives:
    ncount += len(hnative_scores[native])
    

print("Number of specmill sequences found by Hypedsearch", scount, "out of", len(correct_sequences), scount/len(correct_sequences)*100, "%")
print("Number of specmill natives found by Hypedsearch", nacount, "out of", len(specmill_natives), nacount/len(specmill_natives)*100, "%")
print("Number of specmill hybrids found by Hypedsearch", hycount, "out of", len(specmill_hybrids), hycount/len(specmill_hybrids)*100, "%")
print("For Hybrids:\n")
print("Number of Hypedsearch Hybrids", hcount)
print("Number of unique Hypedsearch Hybrids", len(hhybrids))
print("Of these,", hcounth, "were found by SpectrumMill", hcounth/len(hhybrids)*100, "%")
print("Number of Identified Hybrid sequences", hhyb_num)
# print("Number of the top Hybrids that matched with SpectrumMill:", hucounth, hucounth/len(huhybrids)*100, "%")
print("\nFor Natives:\n")
print("Number of Hypedsearch natives", ncount)
print("Number of unique Hypedsearch natives", len(hnatives))
print("Of these,", hcountn, "were found by SpectrumMill", hcountn/len(hnatives)*100, "%")
print("Number of Identified Native sequences", hnat_num)
# print("Number of the top natives that matched with SpectrumMill:", hucountn, hucountn/len(hunaturals)*100, "%")
print("\nData on SpecMill:\n")
print("Number of sequences that SpecMill skipped on the first pass", len(spectra)-len(first_pass_specmill_seqs))
print("Number of sequences that SpecMill skipped on the second pass", len(spectra)-len(correct_sequences))
print("Number of times SpecMill changed it's mind:", changednum)

# Venn diagrams:
#-Specmill naturals vs Hypedsearch naturals
natives_intersection = [value for value in specmill_natives if value in hnatives]
#-Specmill new hybrids vs Hypedsearch top hybrids
hybrids_intersection = [value for value in specmill_hybrids if value in set(hhybrids)]


figure, axes = plt.subplots(1, 2)
v1 = venn2(subsets = (len(specmill_natives), len(hnatives), len(natives_intersection)), set_labels = ('SpectrumMill Natives', 'Hypedsearch Natives'), ax=axes[0])
v2 = venn2(subsets = (len(specmill_hybrids), len(set(hhybrids)), len(hybrids_intersection)), set_labels = ('SpectrumMill new Hybrids', 'Hypedsearch unique Hybrids'), ax=axes[1])
plt.savefig("SMVSHypedsearch")

figure, axis = plt.subplots(1,2)
figure.tight_layout()
#for each hybrid identified by both, print specmill scores in blue and hypedsearch scores in red
max_specmill_intersection_scores = [float(max(native_scores[value])) for value in natives_intersection]
max_hypedsearch_intersection_scores = [float(max(hnative_scores[value])) for value in natives_intersection]
axis[0].scatter(range(0,len(natives_intersection)), max_specmill_intersection_scores, color = 'b')
axis[0].scatter(range(0,len(natives_intersection)), max_hypedsearch_intersection_scores, color = 'r')
axis[0].yaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel("Index")
plt.ylabel("scores")
axis[0].set_title("Specmill and Hypedsearch scores on Natives")

#for each hybrid identified by both, print specmill scores in blue and hypedsearch scores in red
max_specmill_intersection_scores = [float(max(hybrid_scores[value])) for value in hybrids_intersection]
max_hypedsearch_intersection_scores = [float(max(hhybrid_scores[value])) for value in hybrids_intersection]
X_axis = [x for x in range(0,len(hybrids_intersection))]
axis[1].scatter(X_axis, max_specmill_intersection_scores, color = 'b')
axis[1].scatter(X_axis, max_hypedsearch_intersection_scores, color = 'r')
axis[1].yaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel("Index")
plt.ylabel("scores")
axis[1].set_title("Specmill and Hypedsearch scores on Hybrids")
plt.savefig("Score_comparison")

hypedbrids = [x for x in hhybrids]
# figure, axis = plt.subplots()
# hypedsearch_scores, hypedsearch_numbers = get_score_map(hypedbrids, hhybrid_scores)
# axis.bar(hypedsearch_scores, hypedsearch_numbers)
# plt.xlabel("Scores")
# plt.ylabel("Number of hybrids with that score")
# plt.title("Number of hybrids for each score")
# plt.savefig("Number of Hybrids for each score")

figure, axis = plt.subplots()
hypedsearch_scores, hypedsearch_numbers = get_score_map(set(hypedbrids), hhybrid_scores)
axis.bar(hypedsearch_scores, hypedsearch_numbers)
plt.xlabel("Scores")
plt.ylabel("Number of unique hybrids with that score")
plt.title("Number of unique hybrids for each score")
plt.savefig("Number of unique Hybrids for each score")
