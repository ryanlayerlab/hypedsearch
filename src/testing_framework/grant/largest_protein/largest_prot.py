import os
import sys


module_path = os.path.abspath(os.path.join('hypedsearch'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
dataset_path = os.path.join(module_path, "data/database/mouse_database.fasta")

spec_num = -1
max_count = 0
count = 0
max_spec_num = 0
with open(dataset_path, 'r') as d:
    for line in d:
        if line[0] == ">":
            spec_num = spec_num + 1
            if count > max_count:
                max_count = count
                max_spec_num = spec_num
                max_spec_title = line
            count = 0
        else:
            count = count + 1

    print(max_count, max_spec_num)

spec_num = -1
with open(dataset_path, 'r') as d:
    for line in d:
        if line[0] == ">":
            if spec_num == max_spec_num:
                print(line)
                spec_num = spec_num + 1
            else:
                spec_num = spec_num + 1
