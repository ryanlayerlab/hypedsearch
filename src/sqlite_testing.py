from sqlite import database_file
import utils
import gen_spectra
from preprocessing import preprocessing_utils
import os
from testing_framework import testing_utils

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = 0.1
db = database_file(10, False)
# db.check_sizes()

datasets = testing_utils.define_data()
dataset = datasets[0]
input_spectra_path = [os.path.join(dataset[0], 'NOD2_E3.mzML')]
input_spectra, boundaries = preprocessing_utils.load_spectra(input_spectra_path, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)

# query=input("enter query: ")
# db.run_query(query)
# db.query_ion_mass(132.04776143499998, utils.ppm_to_da(132.04776143499998, 10), 'b')
# db.index_ion_mass()
# print('Finished')
input_spectrum = input_spectra[0]
masses = input_spectrum.mz_values

[db.count_ion_mass(target_mass, utils.ppm_to_da(target_mass, 10), 0) for target_mass in masses]