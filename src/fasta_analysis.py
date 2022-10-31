import main

database_file_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/UniProt_mouse.fasta"
count = 0
with open(database_file_path, 'r') as d:
    for line in d:
        if line[0] == '>':
            count = count + 1
        
# db = main.get_database_file(database_file_path)

# protein_list = db.proteins
# total_proteins = len(protein_list)

# total_prot_chars = sum(len(protein_list[key][0][1]) for key in protein_list.keys())
    
print(count)