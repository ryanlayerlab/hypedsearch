import os

module_path = os.path.abspath(os.path.join('hypedsearch'))
filepath = os.path.join(module_path, "data", "database", "mouse_database.fasta")
print(filepath)

count = 0
with open(filepath, 'r') as r:
    for line in r:
        if line[0] == ">":
            count = count + 1
print(count)