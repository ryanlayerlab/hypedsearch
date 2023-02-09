import database
import os
from preprocessing import preprocessing_utils, clustering, merge_search
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import modified_match_masses
import matplotlib.pyplot as plt
from sqlite import database_file
from gen_spectra import get_raw_mass
import matplotlib.pyplot as plt
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25

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
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split('\t')
                correct_sequences.append(split_line[21])
                
    return correct_sequences
def first_pass_truth_set(filepath):
    specmill_seqs = []
    with open(filepath, 'r') as truth_set:
        for line in truth_set:
            split_line = line.split('\t')
            specmill_seqs.append(split_line[7])
    return specmill_seqs
            
def find_specmill_naturals(specmill_path, first_pass_path):
    naturals1, naturals2 = [], []
    with open(specmill_path, 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                if "Hybrid" not in line:
                    split_line = line.split('\t')
                    naturals1.append(split_line[21])
    with open(first_pass_path, 'r') as truth_set:
        for line in truth_set:
            if "Hybrid" not in line:
                split_line = line.split('\t')
                naturals2.append(split_line[7])
    return naturals1, naturals2

specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/AdjustedBMEM_searches.txt'
correct_sequences = truth_set(specmill_path)
first_pass_specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/BMEM_searches.txt'
first_pass_correct_sequences = first_pass_truth_set(first_pass_specmill_path)
first_pass_nat_seqs, second_pass_nat_seqs = find_specmill_naturals(specmill_path, first_pass_specmill_path)

#Things to test:

#1 How many times did spectrumMill find the same thing?
def hypedsearch_top_answers(output_files):
    hhybrids, hnaturals = [], []
    for file in output_files:
        with open(file, 'r') as f:
            prev_id = -1
            for line in f:
                A = line.split('\t')
                if A[0] != prev_id:
                    prev_id = id
                    if A[1] == "Hybrid":
                        hypedsearch_sequence = A[2].replace("-", "")
                        hhybrids.append(hypedsearch_sequence)
                    else:
                        hypedsearch_sequence = A[2]
                        hnaturals.append(hypedsearch_sequence)
    return hnaturals, hhybrids

def hypedsearch_all_answers(output_files):
    hhybrids, hnaturals = [], []
    for file in output_files:
        with open(file, 'r') as f:
            prev_id = -1
            for line in f:
                A = line.split('\t')
                if A[0] != prev_id:
                    prev_id = id
                if A[1] == "Hybrid":
                    hypedsearch_sequence = A[2].replace("-", "")
                    hhybrids.append(hypedsearch_sequence)
                else:
                    hypedsearch_sequence = A[2].replace("-", "")
                    hnaturals.append(hypedsearch_sequence)
    return hnaturals, hhybrids

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
    
open_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/output"
output_files = []
for (root, _, filenames) in os.walk(open_path):
    for fname in filenames:
        output_files.append(os.path.join(root, fname))
hnaturals, hhybrids = hypedsearch_all_answers(output_files)
hunaturals, huhybrids = hypedsearch_top_answers(output_files)

def compare_with_SpecMill(hypedsearch_answer_set, correct_sequences):
    hcount = 0
    for answer in hypedsearch_answer_set:
        if answer in correct_sequences:
            hcount = hcount + 1
    return hcount

def compare_unique_with_SpecMill(hypedsearch_answer_set, correct_sequences):
    hcount = 0
    for answer in hypedsearch_answer_set:
        for specanswer in correct_sequences:
            if test_isomorphic_seqs(answer, specanswer):
                hcount = hcount + 1
    return hcount
    
hcountn = compare_with_SpecMill(hnaturals, correct_sequences)
hcounth = compare_with_SpecMill(hhybrids, correct_sequences)
hucounth = compare_with_SpecMill(huhybrids, correct_sequences)
hucountn = compare_with_SpecMill(hunaturals, correct_sequences)

scount = 0
for answer in correct_sequences:
    if answer in hnaturals or answer in hhybrids:
        scount = scount + 1

print("Number of specmill sequences found by Hypedsearch", scount, "out of", len(correct_sequences), scount/len(correct_sequences)*100, "%")
print("For Hybrids:\n")
print("Number of Hypedsearch Hybrids", len(hhybrids))
print("Of these,", hcounth, "were found by SpectrumMill", len(hhybrids), hcounth/len(hhybrids)*100, "%")
print("Number of Identified Hybrid sequences", len(huhybrids))
print("Number of the top Hybrids that matched with SpectrumMill:", hucounth, hucounth/len(huhybrids)*100, "%")
print("\nFor Naturals:\n")
print("Number of Hypedsearch Naturals", len(hnaturals))
print("Of these,", hcountn, "were found by SpectrumMill", len(hnaturals), hcountn/len(hnaturals)*100, "%")
print("Number of Identified Natural sequences", len(hunaturals))
print("Number of the top Naturals that matched with SpectrumMill:", hucountn, hucountn/len(hunaturals)*100, "%")
print("\nData on SpecMill:\n")
print("Number of sequences that SpecMill skipped on the first pass", len(spectra)-len(first_pass_correct_sequences))
print("Number of sequences that SpecMill skipped on the second pass", len(spectra)-len(correct_sequences))
print("Number of natural sequences in the first pass but not in the second pass", len(first_pass_nat_seqs) - len(second_pass_nat_seqs))