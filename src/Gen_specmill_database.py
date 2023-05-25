import database
from main import get_spectra_files
from preprocessing.preprocessing_utils import load_spectra
import re
import sys

ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 10

dataset_id = 'BMEM_searches'    


prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/UniProt_mouse.fasta'
proteins = database.build(prot_path)

def build_specmill_table(filepath):
    sequences, ids = [], []
    with open(filepath, 'r') as t:
        for line in t:
            A = line.split('\t')
            seq = A[21]
            id = A[14].strip('\"')
            if "BMEM" in id:
                print(A)
            sequences.append(seq)
            ids.append(id)
    return sequences, ids
    
truth_set_pathway = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/' + dataset_id + '.txt'
specmill_sequences, seq_ids = build_specmill_table(truth_set_pathway)


def get_proteins(protein_list, seq_ids, sequences):
    proteins = []
    for i, id in enumerate(seq_ids):
        print("On i=", i)
        found = False
        for protein in protein_list:
            description = protein[0]
            if id in description:
                if sequences[i] in protein[1]:
                    proteins.append(protein)
                    found = True
                    continue
        if 'HYBRID' in id:
            print(id) #these need to be added manually
            continue
            #Ex: HYBRID: mouse ins1C PQVEQLELGGSPGDLQTLAL-DSDKGQQDGFEATTEGPRPQ mouse chgA'
        elif not found:
            print("Not found at", i)

    print(len(seq_ids), len(proteins))
    new_proteins = set(proteins)
    return new_proteins
            
            

all_proteins = get_proteins(proteins.proteins, seq_ids, specmill_sequences)

database_pathway = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/' + dataset_id + '_database.fasta'
with open(database_pathway, 'w') as d:
    for protein in all_proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        seq_len = len(prot_seq)
        curr_write = 0
        prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
        d.write(prot_seq + "\n\n")
