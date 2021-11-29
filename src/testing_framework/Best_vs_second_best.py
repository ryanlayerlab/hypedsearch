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
TOP_X = 50

input_spectra_path = dataset[0]
input_spectra, boundaries = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)

correct_sequences = testing_utils.generate_truth_set(datasets[0])

def calc_diff_array(good_best_array, bad_best_array, good_second_array, bad_second_array):
    good_diff_array = []
    bad_diff_array = []
    for i in range(0, len(good_best_array)):
        target_good = good_best_array[i] - good_second_array[i]
        good_diff_array.append(target_good)
    for i in range(0, len(bad_best_array)):
        target_bad = bad_best_array[i] - bad_second_array[i]
        bad_diff_array.append(target_bad)
    return good_diff_array, bad_diff_array

# def gen_indices(good_diff_array, bad_diff_array):
#     good_indices = []
#     bad_indices = []
#     total_arra
#     return good_indices, bad_indices
    
def plot_best_vs_second_best_for_scoring(good_best_array, bad_best_array, good_second_array, bad_second_array, good_index_array, bad_index_array):
    fig1, ax1 = plt.subplots()
    good_diff_array, bad_diff_array = calc_diff_array(good_best_array, bad_best_array, good_second_array, bad_second_array)
    # good_index_array, bad_index_array = gen_indices(good_diff_array, bad_diff_array)
    ax1.bar(good_index_array, good_diff_array, color = 'g', label = 'Top hit agrees with SpectrumMill')
    ax1.bar(bad_index_array, bad_diff_array, color = 'r', label = 'Top hit does not agree with SpectrumMill')
    plt.title('Diff between best and second best per spectrum')
    plt.xlabel('spectrum num')
    plt.ylabel('diff between top hit and second best')
    plt.legend()
    plt.savefig("Diff between best and second best")

print("Collecting data...")
filepath = os.path.join('..', 'hypedsearch', 'src', 'testing_framework', 'data', 'total_data.txt') #For running locally
# filepath = os.path.join('data', 'total_data.txt') #For running virtually

top_good_hit = []
top_good_prob = []
good_second_array = []
good_second_prob = []
top_bad_hit = []
top_bad_prob = []
bad_second_array = []
bad_second_prob = []
good_indices = []
bad_indices = []
get_ready_for_true_second = False
get_ready_for_false_second = False
with open(filepath, 'r') as d:
    for line in d:
        A = line.rstrip().split('\t')
        spectrum_num = int(A[0])
        score = int(A[1])
        post_prob = float(A[2])
        seq = A[3]
        rank = int(A[4])
        assessment = A[5]
        ion = A[6]

        if get_ready_for_true_second:
            good_second_array.append(score)
            good_second_prob.append(post_prob)
            get_ready_for_true_second = False

        if get_ready_for_false_second:
            bad_second_array.append(score)
            bad_second_prob.append(post_prob)
            get_ready_for_false_second = False

        if assessment == 'True' and rank == 0:
            top_good_hit.append(score)
            top_good_prob.append(post_prob)
            get_ready_for_true_second = True
            good_indices.append(spectrum_num)

        if assessment == 'False' and rank == 0:
            top_bad_hit.append(score)
            top_bad_prob.append(post_prob)
            get_ready_for_false_second = True
            bad_indices.append(spectrum_num)

plot_best_vs_second_best_for_scoring(top_good_hit, top_bad_hit, good_second_array, bad_second_array, good_indices, bad_indices)
# plot_best_vs_second_best_for_scoring(top_good_prob, top_bad_prob, good_second_prob, bad_second_prob, good_indices, bad_indices)