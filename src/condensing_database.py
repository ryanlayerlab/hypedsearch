import database
from main import get_spectra_files
from preprocessing.preprocessing_utils import load_spectra

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 10

dataset_id = 'BMEM_AspN_Fxn4'


prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/UniProt_mouse.fasta'
proteins = database.build(prot_path)


spectra_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/spectra/Lab_Data'
spectra_files = get_spectra_files(spectra_path)
spectra, _ = load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)

def build_specmill_table(filepath):
    sequences, ids = [], []
    with open(filepath, 'r') as t:
        for line in t:
            A = line.split('\t')
            seq = A[7]
            id = A[12][:-1]
            sequences.append(seq)
            ids.append(id)
    return sequences[1:], ids[1:]
    
truth_set_pathway = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/' + dataset_id + '.txt'
specmill_sequences, seq_ids = build_specmill_table(truth_set_pathway)


def get_proteins(protein_list, seq_ids, sequences):
    proteins = []
    for i, id in enumerate(seq_ids):
        for protein in protein_list:
            description = protein[0]
            if id in description: #check for hybrid as well
                if sequences[i] in protein[1]:
                    proteins.append(protein)
                    continue
    
    print(len(seq_ids), len(proteins))
    new_proteins = set(proteins)
    return new_proteins
            
            

all_proteins = get_proteins(proteins.proteins, seq_ids, specmill_sequences)

database_pathway = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/' + dataset_id + '_database.fasta'
with open(database_pathway, 'w') as d:
    for protein in all_proteins:
        d.write(protein[0] + '\n')
        prot_seq = protein[1]
        seq_len = len(prot_seq)
        curr_write = 0
        while seq_len >= 50:
            d.write(prot_seq[curr_write:curr_write+50] + '\n')
            seq_len = seq_len-50
            curr_write = curr_write + 50
        d.write(prot_seq[curr_write:seq_len] + '\n')
        d.write('\n')