import sys, os

module_path = os.path.abspath(os.path.join('..', '..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

import database
import matplotlib.pyplot as plt


# filepath = os.path.abspath(os.path.join('data/database/sample_database.fasta')) #Running locally
filepath = os.path.abspath(os.path.join('../../../data/database/sample_database.fasta')) #Running virtually

def get_subStrings(Str,n):
    substrings = []
    # Pick starting point
    for Len in range(1,n + 1):
         
        # Pick ending point
        for i in range(n - Len + 1):
             
            # Print characters from current
            # starting point to current ending
            # point.
            j = i + Len - 1
 
            for k in range(i,j + 1):
                substrings.append(Str[k] + "")

    return substrings             

def get_kmers(protein_list):
    kmer_dict = dict()
    for prot in protein_list:
        seq = protein_list[prot][0].sequence
        substrings = get_subStrings(seq, len(seq))
        for str in substrings:
            if len(str) not in kmer_dict:
                kmer_dict[len(str)] = 1
            else:
                kmer_dict[len(str)] = kmer_dict[len(str)] + 1
    return kmer_dict

def plot_occurance(output_dir, kmer_dict):
    
    plt.plot(kmer_dict.keys(), [kmer_dict[x] for x in kmer_dict.keys()])
    plt.title('Occurance by kmer length')
    plt.xlabel('Len of kmer')
    plt.ylabel('Occurance')
    plt.savefig(output_dir)

db = database.build(filepath)
print('Building kmer dictionary')
kmer_dict = get_kmers(db.proteins)
print('Done')
print('Plotting')
output_dir = os.path.abspath(os.path.join('../plots/Occurance.txt'))
plot_occurance(output_dir, kmer_dict)
print('Done')