import database
from preprocessing import preprocessing_utils, clustering, merge_search
from main import get_spectra_files
from utils import ppm_to_da
from preprocessing.merge_search import modified_match_masses
import matplotlib.pyplot as plt
from sqlite import database_file
from gen_spectra import get_raw_mass
import matplotlib.pyplot as plt
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25

prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta'
proteins = database.build(prot_path)

spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/NOD2_E3'
spectra_files = get_spectra_files(spectra_path)
spectra, boundaries = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

def truth_set(filepath):
    correct_sequences = []
    with open(filepath, 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split(';')
                correct_sequences.append(split_line[9])
                
    return correct_sequences

specmill_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/NOD2_E3_results.ssv'
correct_sequences = truth_set(specmill_path)

#Optimizations I can make

#Idea 1: Only consider unique hits
#Q: How to determine uniqueness without the string - IDK
#Q: How much does the number of b_clusters and y_clusters change when we only consider unique hits.
def find_sequence(pid, start, end, protein_list):
    prot_seq = protein_list[pid][1]
    return prot_seq[start:end]

def get_unique_clusters(clusters, protein_list):
    unique_clusters = set()
    for A in clusters:
        pid = int(A[1])
        start = int(A[3])
        end = int(A[4])
        seq = find_sequence(pid, start, end, protein_list)
        unique_clusters.add(seq)
    return unique_clusters

def create_hits(spec_num,spectrum,matched_masses_b,matched_masses_y):
    b_hits, y_hits = [], []
    for mz in spectrum.mz_values:
        if mz in matched_masses_b:
            for tuple in matched_masses_b[mz]:
                tup = (spec_num, mz, tuple)
                b_hits.append(tup)
        if mz in matched_masses_y:
            for tuple in matched_masses_y[mz]:
                tup = (spec_num, mz, tuple)
                y_hits.append(tup)
    return b_hits, y_hits

percentage_list = []
for spectrum in spectra:
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, False)
    b_hits,y_hits = create_hits(spectrum.num,spectrum,matched_masses_b,matched_masses_y)
    clusters = clustering.create_clusters('b', b_hits, y_hits)
    cluster_size = len(clusters)
    unique_clusters = get_unique_clusters(clusters, proteins.proteins)
    unique_size = len(unique_clusters)
    print(spectrum.num, "Unique clusters are", unique_size/cluster_size * 100, "percent the normal size", unique_size, cluster_size)
    percentage_list.append(unique_size/cluster_size * 100)

index_list = []
[index_list.append(i) for i in range(0,len(spectra))]
plt.scatter(index_list, percentage_list)
plt.xlabel("Spectrum Index")
plt.ylabel("Percentage of the normal dataset size")
plt.title
plt.savefig("Percent reduction with unique kmers")

print("On average, unique kmers are", sum(percentage_list)/len(percentage_list), 'percent the size of the normal dataset')

#Idea 2: Reduce database size (by a minor factor) by forcing the charge on the peaks that we can

#Idea 3: Try and figure out if we can deduce a hybrid from the spectrum and it's match with the precursor

#Idea 4: Dynamic programming to fid sequences that match the precursor mass