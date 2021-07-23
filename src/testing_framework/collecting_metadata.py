import os
import matplotlib.pyplot as plt
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

# refresh metadata
with open('metadata.txt', 'w') as m:
    m.close()

datasets = testing_utils.define_data()

correct_sequences = testing_utils.generate_truth_set(datasets[0])
# BALB3_set = testing_utils.generate_truth_set(BALB3_data)
# [correct_sequences.append(x) for x in BALB3_set]

ppm_tolerance = 20

for dataset in datasets:
    input_spectra_path = dataset[0]
    spectra = []
    input_list = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)
    [spectra.append(x) for x in input_list]

# refresh metadata
with open('metadata.txt', 'w') as m:
    m.close()

print('Collecting metadata...')
testing_utils.collect_metadata()
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
ListOfLen = []

for i, spectrum in enumerate(spectra):
    max_abundance = 0
    found = False
    initial_hits = []
    input_spectrum = spectrum[0]
    tot_measured_spec_length = tot_measured_spec_length + len(input_spectrum)
    input_abundance_set = spectrum[1]
    precursor_mass = spectrum[5]
    precursor_charge = spectrum[6]
    ideal_spectrum = gen_spectra.gen_spectrum(correct_sequences[i])
    tot_ideal_spec_length = tot_ideal_spec_length + len(ideal_spectrum['spectrum'])
    # Checking input_spectrum for hits
    for a, j in enumerate(input_spectrum):
        #Finding max abundance
        if (input_abundance_set[a] > max_abundance):
            max_abundance = input_abundance_set[a]
            max_abundance_location = a

        for k in ideal_spectrum['spectrum']:
            tolerance = utils.ppm_to_da(j, ppm_tolerance)
            if testing_utils.isintolerance(j, k, tolerance):
                initial_hits.append(j)
                all_hits.append(j)
                hit_abundances.append(input_abundance_set[a])
                found = True

        if found == False:
            all_misses.append(j)
            miss_abundances.append(input_abundance_set[a])

    misleading = False
    for l in ideal_spectrum['spectrum']:
        tolerance = utils.ppm_to_da(l, ppm_tolerance)
        if testing_utils.isintolerance(input_spectrum[max_abundance_location], l, tolerance):
            misleading = True

    if (misleading == True):
        misleading_abundances.append(i)

    # Checking precursor mass
    tolerance = utils.ppm_to_da(precursor_mass, ppm_tolerance)
    if (not testing_utils.isintolerance(precursor_mass, gen_spectra.get_precursor(correct_sequences[i], precursor_charge), tolerance)):
        count = count + 1

    with open('metadata.txt', 'a') as m:
        m.write(str(i) + ': ')
        m.write(str(initial_hits) + '\n')
        ListOfLen.append(len(initial_hits))
        total_hits = total_hits + len(initial_hits)

avg_hits = total_hits / i
# Collecting average abundance
avg_hit_abundance = testing_utils.get_average_from_set(hit_abundances)
avg_miss_abundance = testing_utils.get_average_from_set(miss_abundances)
# Collecting total length
total_length = testing_utils.get_total_length(correct_sequences)

#Writing data
with open('metadata.txt', 'a') as m:
    m.write('total number of hits: ' + str(len(all_hits)) + ' out of ' + str(tot_ideal_spec_length) + ' potential hits ' + '(' + str(round((len(all_hits) / tot_ideal_spec_length) * 100)) + '%)' + '\n')
    m.write('total number of misses: ' + str(len(all_misses)) + ' out of ' + str(tot_ideal_spec_length - len(all_hits)) + ' potential misses ' + '(' + str(round((len(all_misses) / tot_ideal_spec_length) * 100)) + '%)' + '\n')
    m.write(str(round((total_hits/tot_measured_spec_length) * 100)) + '% of m/z values in input spectra correlated to a correct hit \n')
    m.write(str(round((len(all_misses)/tot_measured_spec_length) * 100)) + '% of m/z values in input spectra correlated to a miss \n')
    m.write('average number of hits per spectrum: ' + str(avg_hits) + '\n')
    m.write('average length of correct sequences: ' + str(total_length / i) + '\n')
    m.write('Precursor was not found in data for ' + str(count) + '/' + str(i + 1) + ' spectra\n')
    m.write('Average abundance of a hit: ' + str(avg_hit_abundance) + '\n')
    m.write('Average abundance of a miss: ' + str(avg_miss_abundance) + '\n')
    m.write(str(len(misleading_abundances)) + ' times it the largest abundance was misleading, this happened at ' + str(misleading_abundances) + '\n')
    m.write('max abundance of a hit: ' + str(max(hit_abundances)) + ', min abundance of a hit: ' + str(min(hit_abundances)) + '\n')
    m.write('max abundance of a miss: ' + str(max(miss_abundances)) + ', min abundance of a miss: ' + str(min(miss_abundances)) + '\n')

print('Done')

arr = np.array(all_hits)

plt.scatter(all_hits, hit_abundances, s=None, c='green')
plt.scatter(all_misses, miss_abundances, s=None, c='red')
plt.show(block=True)