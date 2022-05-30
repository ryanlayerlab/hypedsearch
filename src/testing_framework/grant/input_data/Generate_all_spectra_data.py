import os
import sys

# module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../../..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from preprocessing import preprocessing_utils
 
ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
peak_filter = 25
relative_abundance_filter = 0.1
size_lim = 10
cap_size = 20
get = True
no_k = True

database_path = os.path.abspath(os.path.join('../../../../data/spectra/NOD2_E3'))
input_spectra_path = [os.path.join(database_path, 'NOD2_E3.mzML')]
input_spectra, boundaries = preprocessing_utils.load_spectra(input_spectra_path, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)

for i, spectrum in enumerate(input_spectra):
    with open("Spec_" + str(i) + "_input_data.txt", "w") as s:
        s.write(str(spectrum.precursor_mass)+"\t"+str(spectrum.precursor_charge)+"\n")
        for i, x in enumerate(spectrum.mz_values):
            abundance = spectrum.abundance[i]
            s.write(str(x) + "\t" + str(abundance) + "\n")