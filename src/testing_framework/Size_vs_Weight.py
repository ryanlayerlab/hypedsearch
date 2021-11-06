import matplotlib.pyplot as plt
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

import testing_utils
import gen_spectra

datasets = testing_utils.define_data()

dataset = datasets[0]
ppm_tolerance = 20

input_spectra_path = dataset[0]
input_spectra, boundaries = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)

correct_sequences = testing_utils.generate_truth_set(datasets[0])


def plot_size_vs_weight(good_size_array, bad_size_array, good_weight_array, bad_weight_array):
    fig1, ax1 = plt.subplots()
    ax1.scatter(bad_size_array, bad_weight_array, color = 'r', label = 'No hit agrees with SpectrumMill in the top 50')
    ax1.scatter(good_size_array, good_weight_array, color = 'g', label = 'Hit agrees with SpectrumMill in top 50')
    plt.title('Size vs Weight')
    plt.xlabel('Size')
    plt.ylabel('Weight')
    plt.legend()
    plt.savefig("Size vs Weight")

def see_if_good_hit_in_top_50(b_cluster_array, y_cluster_array, SpecMill_sequence):
    for cluster in b_cluster_array:
        seq = cluster[3]
        result, _ = testing_utils.is_good_hit(seq, 'b', SpecMill_sequence)
        return result


print("Collecting data...")
filepath = os.path.join('..', 'hypedsearch', 'src', 'testing_framework', 'data', 'total_data.txt') #For running locally
# filepath = os.path.join('data', 'total_data.txt') #For running virtually
b_cluster_array = []
y_cluster_array = []
good_size_array = []
bad_size_array = []
good_weight_array = []
bad_weight_array = []
prev_spectrum_num = 0
with open(filepath, 'r') as d:
    for i, line in enumerate(d):
        A = line.rstrip().split('\t')
        spectrum_num = int(A[0])
        score = int(A[1])
        post_prob = float(A[2])
        seq = A[3]
        rank = int(A[4])
        assessment = A[5]
        ion = A[6]

        if ion == 'b':
            b_cluster_array.append((spectrum_num, score, post_prob, seq, rank, assessment, ion))
        else:
            y_cluster_array.append((spectrum_num, score, post_prob, seq, rank, assessment, ion))
        if spectrum_num != prev_spectrum_num:
            if see_if_good_hit_in_top_50(b_cluster_array, y_cluster_array, correct_sequences[spectrum_num]):
                good_size_array.append(len(correct_sequences[spectrum_num]))
                good_weight_array.append(gen_spectra.get_precursor(correct_sequences[spectrum_num]))
            else:
                bad_size_array.append(len(correct_sequences[spectrum_num]))
                bad_weight_array.append(gen_spectra.get_precursor(correct_sequences[spectrum_num]))
            cluster_array = []
            prev_spectrum_num = spectrum_num

plot_size_vs_weight(good_size_array, bad_size_array, good_weight_array, bad_weight_array)