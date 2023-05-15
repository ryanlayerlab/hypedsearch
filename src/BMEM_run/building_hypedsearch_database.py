import os
import yaml
# Comet seems like the best scored one goes on top but they don't use integer scoring
# Want code that takes in the raw_Comet_results and builds a suitable database for hypedsearch to use

def get_files(folderpath):
    file_list = []
    for (root, _, filenames) in os.walk(folderpath):
        for fname in filenames:
            file_list.append(os.path.join(root, fname))
    return file_list

def get_raw_comet_results(folderpath):
    print(folderpath)
    files = sorted(get_files(folderpath))
    results_dict = dict()
    for file in files:
        results_dict[file] = []
    
    return results_dict

folderpath = os.path.join(os.path.abspath(__file__), "..", "raw_Comet_results")
results = get_raw_comet_results(folderpath)
