import os, re
from itertools import product
import numpy as np

HYBRID_ALIGNMENT_PATTERN = re.compile(r'[-\(\)]')

def string_to_bool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def file_exists(file_name: str) -> bool:
    return os.path.isfile(file_name)

def make_valid_dir_string(dir_path: str) -> str:
    return f'{dir_path}{os.path.sep}' if not dir_path.endswith(os.path.sep) else dir_path

def make_dir(dir_path: str) -> bool:
    dir_path = make_valid_dir_string(dir_path)
    try:
        if not os.path.exists(dir_path): 
            os.makedirs(dir_path)
        return True 
    except:
        return False

def make_valid_json_file(file_name: str) -> str:
    if not file_name.endswith('.json'): file_name += '.json'
    return file_name

def make_valid_fasta_file(file_name: str) -> str:
    if not file_name.endswith('.fasta'): file_name += '.fasta'
    return file_name

def is_json(file: str) -> bool:
    return file.endswith('.json')

def is_fasta(file: str) -> bool:
    return file.endswith('.fasta')

def is_file(file: str) -> bool:
    return os.path.isfile(file)

def all_perms_of_s(s: str, keyletters: str) -> list:
    seq = list(s)
    perms = []
    keyletter_indices = [ i for i, c in enumerate(seq) if c in keyletters ]
    for t in product(keyletters, repeat=len(keyletter_indices)):
        for i, c in zip(keyletter_indices, t):
            seq[i] = c
        perms.append(''.join(seq))
    return perms

def ppm_to_da(mass: float, ppm_tolerance: float) -> float:
    return abs((ppm_tolerance / 1000000)*mass)

def make_sparse_array(spectrum: list, width: float, value=50) -> np.ndarray:
    list_size = int(max(spectrum)//width) + 1
    sparse = np.zeros(list_size)
    for m in spectrum:
        sparse[int(m // width)] = value
    return sparse

def to_percent(index, total):
    return int(100 * (index)/total)

def hashable_boundaries(boundaries: list) -> str:
    if (len(boundaries) == 2):
        return '-'.join([str(x) for x in boundaries])
    else:
        return None



