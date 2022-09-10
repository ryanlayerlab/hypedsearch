# scale linearly
# based on not protein list length
# but based on the sum of lengths of all proteins


# ABC
# AB
# ABCDE


# output -> number of proteins and total sum of lengths of 

# ABJKASHBDKUASHDKJASHDJASGDKASYDJGASDASGDKASUYDBHASJDMSAD

# prot 274 -> 75 characters, 11K records -> ratio X
# sample databse -> N chars, NX records
# mouse databse -> M chars, MX records

import main

database_file_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/sample_database.fasta"
db = main.get_database_file(database_file_path)

protein_list = db.proteins
total_proteins = len(protein_list)
# total_prot_chars = 0
# for key in protein_list.keys():
#     protein = protein_list[key]
#     sequence = protein[0][1]
#     total_prot_chars = total_prot_chars + len(sequence)

total_prot_chars = sum(len(protein_list[key][0][1]) for key in protein_list.keys())
    
print(total_proteins,total_prot_chars)