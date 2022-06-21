import os
import sys

module_path = os.path.abspath(os.path.join('hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
# module_path = os.path.abspath(os.path.join('../../..'))
# if module_path not in sys.path:
#     sys.path.append(module_path) 

natural_spectra, hybrid_spectra = [], []
data_filepath = os.path.join(module_path, '../data/output/summary.txt')
prev_spec_num = -1
with open(data_filepath, 'r') as r:
    for i, line in enumerate(r):
        A = line.rstrip().split('\t')
        if A[0] != prev_spec_num:
            prev_spec_num = A[0]
            if A[1] == 'Natural':
                natural_spectra.append(prev_spec_num)
            else:
                hybrid_spectra.append(prev_spec_num)

with open('spec_with_top_natural.txt', 'w') as w:
    [w.write(str(x) + "\n") for x in natural_spectra]
with open('spec_with_top_hybrid.txt', 'w') as w:
    [w.write(str(x) + "\n") for x in hybrid_spectra]