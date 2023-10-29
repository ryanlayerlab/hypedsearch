import os
import matplotlib.pyplot as plt
from statistics import mean
import numpy as np
from collections import namedtuple

import sys
module_path = os.path.abspath(os.path.join('..', 'hypedsearch'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)
import testing_utils
import utils, gen_spectra


#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20

# refresh metadata
with open('metadata.txt', 'w') as m:
    m.close()

datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)

correct_sequences = testing_utils.generate_truth_set(datasets[0])
correct_sequence = correct_sequences[0]

print('Collecting metadata...')
avg_hits = 0
total_hits = 0
tot_ideal_spec_length = 0
tot_measured_spec_length = 0
count = 0

hit_abundances = []
miss_abundances = []
misleading_abundances = []
all_hits = []
all_misses = []

testing_utils.collect_metadata(input_spectra, correct_sequences, ppm_tolerance, all_hits, all_misses, hit_abundances, miss_abundances, misleading_abundances)

#Writing data
with open('metadata.txt', 'a') as m:
    m.write('total number of hits: ' + str(len(all_hits)) + ' out of ' + str(tot_ideal_spec_length) + ' potential hits ' + '(' + str(round((len(all_hits) / tot_ideal_spec_length) * 100)) + '%)' + '\n')
    m.write('total number of misses: ' + str(len(all_misses)) + ' out of ' + str(tot_ideal_spec_length - len(all_hits)) + ' potential misses ' + '(' + str(round((len(all_misses) / tot_ideal_spec_length) * 100)) + '%)' + '\n')
    m.write(str(round((total_hits/tot_measured_spec_length) * 100)) + '% of m/z values in input spectra correlated to a correct hit \n')
    m.write(str(round((len(all_misses)/tot_measured_spec_length) * 100)) + '% of m/z values in input spectra correlated to a miss \n')
    m.write('average number of hits per spectrum: ' + str(avg_hits) + '\n') #Replace with distribution
    m.write('average length of correct sequences: ' + str(total_length / i) + '\n') #Replace with distribution
    m.write('Precursor was not found in data for ' + str(count) + '/' + str(i + 1) + ' spectra\n')
    m.write('Average abundance of a hit: ' + str(mean(hit_abundances)) + '\n') #Replace with distribution
    m.write('Average abundance of a miss: ' + str(mean(miss_abundances)) + '\n') #Replace with distribution
    m.write(str(len(misleading_abundances)) + ' times it the largest abundance was misleading, this happened at ' + str(misleading_abundances) + '\n')
    m.write('max abundance of a hit: ' + str(max(hit_abundances)) + ', min abundance of a hit: ' + str(min(hit_abundances)) + '\n')
    m.write('max abundance of a miss: ' + str(max(miss_abundances)) + ', min abundance of a miss: ' + str(min(miss_abundances)) + '\n')

print('Done')