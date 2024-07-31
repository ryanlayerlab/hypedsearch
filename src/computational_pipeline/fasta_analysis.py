import main

database_file_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/UniProt_mouse.fasta"
count = 0
with open(database_file_path, 'r') as d:
    for line in d:
        if line[0] == '>':
            count = count + 1
print(count)