import os, re
from itertools import product
import numpy as np

HYBRID_ALIGNMENT_PATTERN = re.compile(r'[-\(\)]')

def string_to_bool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def boolean_string(s): # not sure what the purpose of this func is. what is it supposed to return exactly?
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'False'

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

def make_valid_text_file(file_name: str) -> str:
    if not file_name.endswith('.txt'): file_name += '.txt'
    return file_name

def make_valid_json_file(file_name: str) -> str:
    if not file_name.endswith('.json'): file_name += '.json'
    return file_name

def make_valid_csv_file(file_name: str) -> str:
    if not file_name.endswith('.csv'): file_name += '.csv'
    return file_name

def make_valid_fasta_file(file_name: str) -> str:
    if not file_name.endswith('.fasta'): file_name += '.fasta'
    return file_name

def is_json(file: str) -> bool:
    return file.endswith('.json')

def is_fasta(file: str) -> bool:
    return file.endswith('.fasta')

def is_dir(dir_path: str) -> bool:
    return os.path.isdir(dir_path)

def is_file(file: str) -> bool:
    return os.path.isfile(file)

def all_perms_of_s(s: str, keyletters: str) -> list:
    seq = list(s)
    perms = []
    indices = [ i for i, c in enumerate(seq) if c in keyletters ]
    for t in product(keyletters, repeat=len(indices)):
        for i, c in zip(indices, t):
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

def cosine_similarity(a: list, b: list) -> float:
    if len(a) > len(b):
        b = list(b) + list(np.zeros(len(a) - len(b), dtype = int))
    elif len(b) > len(a):
        a = list(a) + list(np.zeros(len(b) - len(a), dtype = int))

    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def split_hybrid(sequence: str):
    if '-' in sequence:
        return tuple(sequence.split('-', 2)[:2])
    else:
        left = sequence.split(')', 1)[0].replace('(', '')
        right = sequence.split('(', 2)[1].replace(')', '')
        return (left, right)

def find_dir(filename, location):
    filepath = os.path.join(location, filename)
    return os.path.exists(filepath)

def DEV_contains_truth_parts(truth_seq: str, hybrid: bool, b_seqs: list, y_seqs: list) -> bool:
    replacer_small = lambda x: re.sub('I|L', 'B', x)
    def replacer_big(x):
        return re.sub('-|(|)', '', replacer_small(x))
    
    has_left = has_right = False
    left_half = right_half = ''
    if hybrid:
        if '-' in truth_seq:
            left_half, right_half = truth_seq.split('-', 2)[:2]

        elif '(' in truth_seq and ')' in truth_seq:
            left_half = truth_seq.split(')')[0].replace('(', '')
            right_half = truth_seq.split('(')[1].replace(')', '')

        else: 
            left_half = truth_seq[:2]
            right_half = truth_seq[-2:]

        truth_seq = replacer_big(truth_seq)

        b_seqs = [replacer_small(x) for x in b_seqs]
        y_seqs = [replacer_small(x) for x in y_seqs]
    
    filtered_b_seqs = list(filter(lambda x: len(x) > 1), b_seqs)
    filtered_y_seqs = list(filter(lambda x: len(x) > 1), y_seqs)
    has_left = any(truth_seq.startswith(x) for x in filtered_b_seqs) or \
               any(x.startswith(truth_seq) for x in filtered_b_seqs)

    has_right = any(truth_seq.endswith(x) for x in filtered_y_seqs) or \
                any(x.endswith(truth_seq) for x in filtered_y_seqs)

    if hybrid:
        if not has_left:
            left_half = replacer_small(left_half)
            has_left = any(x.startswith(left_half) for x in b_seqs)

        if not has_right:
            right_half = replacer_small(right_half)
            has_right = any(x.endswith(right_half) for x in y_seqs)

        return has_left and has_right

    return has_left or has_right

def CICD_test():
    return 1

def DEV_contains_truth_exact(truth_seq: str, hybrid: bool, seqs: list) -> bool:
    def replacer(x):
        temp = re.sub('I|L', 'B', x)
        return re.sub('-|(|)', '', temp)
    
    if hybrid:
        truth_seq = replacer(truth_seq)
        seqs = [replacer(x) for x in seqs]

    return truth_seq in seqs



