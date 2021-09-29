import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)


import database
import testing_utils
import operator

from objects import Spectrum

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20



datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)

correct_sequences = testing_utils.generate_truth_set(datasets[0])

path = dataset[2]
db = database.build(path)






matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)
print('Finished matching masses')





print('Collecting data...')
with open('data.txt', 'w') as d:
    d.write('')

for spectrum_num, input_spectrum in enumerate(input_spectra):
    
    
    
    
    
    
    correct_sequence = correct_sequences[spectrum_num]

    # input_spectrum = input_spectra[spectrum_num]

    correct_hits = []
    #Remember to add in abundance if it is helpful
    b_hits, y_hits, b_set, y_set, misses = testing_utils.find_hits(mz_mapping, boundaries, input_spectrum, spectrum_num, matched_masses_b, matched_masses_y)
    testing_utils.append_correct_hits(correct_hits, correct_sequence, input_spectrum, ppm_tolerance)





    #Writing b and y hits
    testing_utils.write_hits(b_hits, y_hits)

    ion = 'b'
    testing_utils.create_clusters(ion)
    b_sorted_clusters = testing_utils.sort_clusters_by_post_prob(ion, mz_mapping, boundaries, matched_masses_b, matched_masses_y)
    testing_utils.write_b_sorted_cluster(b_sorted_clusters)

    ion = 'y'
    testing_utils.create_clusters(ion)
    y_sorted_clusters = testing_utils.sort_clusters_by_post_prob(ion, mz_mapping, boundaries, matched_masses_b, matched_masses_y)
    testing_utils.write_y_sorted_cluster(y_sorted_clusters)

    b_sorted_clusters = sorted(b_sorted_clusters, key=operator.attrgetter('score', 'post_prob', 'pid', 'prior'), reverse = True)
    y_sorted_clusters = sorted(y_sorted_clusters, key=operator.attrgetter('score', 'post_prob', 'pid', 'prior'), reverse = True)

    for i, cluster in enumerate(b_sorted_clusters):
        score = cluster.score
        post_prob = cluster.post_prob
        seq = cluster.seq
        cluster_num = i
        ion = 'b'
        assessment, _ = testing_utils.is_good_hit(cluster.seq, ion, correct_sequence)


        with open("data.txt", 'a') as d:
            d.write(str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')

    for i, cluster in enumerate(y_sorted_clusters):
        score = cluster.score
        post_prob = cluster.post_prob
        seq = cluster.seq
        cluster_num = i
        ion = 'y'
        assessment, _ = testing_utils.is_good_hit(cluster.seq, ion, correct_sequence)
    
        with open("data.txt", 'a') as d:
            d.write(str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')
print('Done')
