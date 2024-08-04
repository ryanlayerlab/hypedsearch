import database
from preprocessing import preprocessing_utils
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import get_all_matched_rows
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

prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(prot_path)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/NOD2_E3'
spectra_files = get_spectra_files(spectra_path)
spectra, boundaries = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

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

def get_raw_prec(precursor, prec_charge):
    # seqmass + water 
    raw_prec = (precursor * prec_charge) - (prec_charge * PROTON_MASS) - WATER_MASS
    return raw_prec

def make_data_structures(total_weight):
    total_string = str(total_weight)
    for i,char in enumerate(total_string):
        if char == '.':
            to_mult = 10**(len(total_string[i:])+1)
    total_string = total_string.replace('.', '')
    converted_weight = int(total_string)
    val = []
    [val.append(1) for x in range(0,len(AMINO_ACIDS))]
    wt = []
    for aa in AMINO_ACIDS:
        wt.append(AMINO_ACIDS[aa]*to_mult)
    
    return val, wt, converted_weight, to_mult 
 
def knapsack(wt, val, W, n, t):
 
    # base conditions
    if n == 0 or W == 0:
        return 0
    if t[n][W] != -1:
        return t[n][W]
 
    # choice diagram code
    if wt[n-1] <= W:
        t[n][W] = max(
            val[n-1] + knapsack(
                wt, val, W-wt[n-1], n-1),
            knapsack(wt, val, W, n-1))
        return t[n][W]
    elif wt[n-1] > W:
        t[n][W] = knapsack(wt, val, W, n-1)
        return t[n][W]
 
 # We can encode all the amino acid weights as prime numbers so that the product uniquely factors back into it's primes

def dissect_hit(result, to_mult, weights):
    result = result/to_mult
    final_seq = ''
    #now that the result is back in decimals, we want to find out the sequence that multiplies back to this
    for aa in AMINO_ACIDS:
        while result % AMINO_ACIDS[aa] == 0:
            final_seq.append(aa)
    return final_seq

def initialize_table(values):
    for value in values:
    
for spectrum in spectra:
    raw_prec_mass = get_raw_prec(spectrum.precursor_mass, spectrum.precursor_charge)    
    W = raw_prec_mass + ppm_to_da(spectrum.precursor_mass, spectrum.precursor_charge)
    vals, weights, conv_W, to_mult = make_data_structures(W)
    print(str(conv_W))
    n = len(vals)
    # We initialize the matrix with -1 at first.
    t = [[-1 for i in range(conv_W + 1)] for j in range(n + 1)]
    t = initialize_table()

    aligned_val = knapsack(weights, vals, conv_W, n, t)
    aligned_seq = dissect_hit(aligned_val, to_mult)
    print(aligned_seq)