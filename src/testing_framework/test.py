import os
import sys
from typing import Tuple
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

import gen_spectra
import testing_utils
from scoring.mass_comparisons import optimized_compare_masses

datasets = testing_utils.define_data()
ppm_tolerance = 20

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


hit = "KNDIPKDK"


# 5	0.8110311825930346	2	152	KNDIPKDK	939.4456176757812	(78, 78, 'K', '129.1007537841797')	(78, 79, 'KN', '243.13357543945312')
# 	(78, 80, 'KND', '358.1739196777344')	(78, 84, 'KNDIPKD', '811.4014892578125')	(78, 85, 'KNDIPKDK', '939.4456176757812')
input_mz = [70.06548309326172, 84.04374694824219, 84.0806655883789, 86.09642028808594, 101.07074737548828, 102.05448913574219, 129.1007537841797, 132.10121154785156, 225.12210083007812, 226.1183624267578, 235.1075439453125, 243.13357543945312, 261.1429443359375, 341.1453857421875, 349.1622619628906, 358.1739196777344, 412.21746826171875, 440.216064453125, 569.2589721679688, 586.2818603515625, 697.3143310546875, 698.3234252929688, 810.3919677734375, 811.4014892578125, 939.4456176757812]
b_spec = gen_spectra.gen_spectrum(hit, charge=None, ion='b')
addscore = optimized_compare_masses(mz_mapping, boundaries, input_mz, b_spec, 20, False)
print(addscore)