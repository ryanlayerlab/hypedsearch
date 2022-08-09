import sys
import os

from sqlalchemy import false

module_path = os.path.abspath(os.path.join('hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
# module_path = os.path.abspath(os.path.join('../../..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)

from gen_spectra import get_precursor
from testing_framework import testing_utils
from preprocessing import preprocessing_utils
import database

ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
peak_filter = 25
relative_abundance_filter = 0.1

def get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter):
    datasets = testing_utils.define_data()
    dataset = datasets[0]
    input_spectra_path = [os.path.join(dataset[0], 'NOD2_E3.mzML')]
    input_spectra, boundaries = preprocessing_utils.load_spectra(input_spectra_path, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
    correct_sequences = testing_utils.generate_truth_set(datasets[0])
    path = dataset[2]
    db = database.build(path)

    return input_spectra, boundaries, correct_sequences, db

input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)

input_precs, input_charges = dict(), dict()
for i, spectrum in enumerate(input_spectra):
    input_precs[i] = spectrum.precursor_mass
    input_charges[i] = spectrum.precursor_charge

hypedsearch_results, specMill_results, hypedsearch_closer, specmill_closer, eq_dist = [],[],[],[],[]
with open("hypedsearch/src/testing_framework/grant/SpectrumMill_precursor_dist/hs_v_sm.txt", "r") as h:
    for line in h:
        A = line.rstrip().split(' ')
        spec_num = int(A[0])
        hypedsearch_result = A[2]
        specmill_result = A[4]
        evaluation = A[1]
        agreement = A[5]
        
        if evaluation == 'Hybrid':
            hypedsearch_result = hypedsearch_result.replace("-", "")
        
        hyped_dist = abs(input_precs[spec_num] - get_precursor(hypedsearch_result, input_charges[spec_num]))
        specmill_dist = abs(input_precs[spec_num] - get_precursor(specmill_result, input_charges[spec_num]))
        
        
        if agreement == 'False':
            if hyped_dist < specmill_dist:
                hypedsearch_closer.append(line)
            elif specmill_dist < hyped_dist:
                specmill_closer.append(line)
            else:
                eq_dist.append(line)
                
with open("hypedsearch/src/testing_framework/grant/SpectrumMill_precursor_dist/hypedsearch_closer.txt", 'w') as c:
    [c.write(x) for x in hypedsearch_closer]

with open("hypedsearch/src/testing_framework/grant/SpectrumMill_precursor_dist/specmill_closer.txt", 'w') as c:
    [c.write(x) for x in specmill_closer]
    
with open("hypedsearch/src/testing_framework/grant/SpectrumMill_precursor_dist/eq_dist.txt", 'w') as c:
    [c.write(x) for x in eq_dist]