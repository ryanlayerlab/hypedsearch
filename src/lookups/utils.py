import os, re
from itertools import product
import numpy as np

HYBRID_ALIGNMENT_PATTERN = re.compile(r'[-\(\)]')

def string_to_bool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'False'

def file_exists(file_name: str) -> bool:
    return os.path.isfile(file_name)

def make_valid_dir_string(dir_path: str) -> str:
    return dir_path + os.path.sep if os.path.sep != dir_path[-1] else dir_path

def make_dir(dir_path: str) -> bool:
    dir_path = make_valid_dir_string(dir_path)
    try:
        if not os.path.exists(dir_path): 
            os.makedirs(dir_path)
        return True 
    except:
        return False

def make_valid_text_file(file_name: str) -> str:
    file_name = file_name + '.txt' if '.txt' not in file_name else file_name
    return file_name

def make_valid_json_file(file_name: str) -> str:
    file_name = file_name + '.json' if '.json' not in file_name else file_name
    return file_name

def make_valid_csv_file(file_name: str) -> str:
    file_name = file_name + '.csv' if '.csv' not in file_name else file_name
    return file_name

def make_valid_fasta_file(file_name: str) -> str:
    file_name = file_name + '.fasta' if '.fasta' not in file_name else file_name
    return file_name

def is_json(file: str) -> bool:
    return True if '.json' in file else False

def is_fasta(file: str) -> bool:
    return True if '.fasta' in file else False

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
        b = list(b) + [0 for _ in range(len(a) - len(b))]
        
    if len(b) > len(a):
        a = list(a) + [0 for _ in range(len(b) - len(a))]

    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def split_hybrid(sequence: str):
    if '-' in sequence:
        return (sequence.split('-')[0], sequence.split('-')[1])
    
    else:
        left = sequence.split(')')[0].replace('(', '')
        right = sequence.split('(')[1].replace(')', '')
        return (left, right)

def find_dir(filename, location):
    filepath = os.path.join(location, filename)
    if os.path.exists(filepath):
        return True
    else:
        return False

def DEV_contains_truth_parts(truth_seq: str, hybrid: bool, b_seqs: list, y_seqs: list) -> bool:
    has_left = False
    has_right = False
    left_half = ''
    right_half = ''
    if hybrid:
        if '-' in truth_seq:
            left_half = truth_seq.split('-')[0]
            right_half = truth_seq.split('-')[1]

        elif '(' in truth_seq and ')' in truth_seq:
            left_half = truth_seq.split(')')[0].replace('(', '')
            right_half = truth_seq.split('(')[1].replace(')', '')

        else: 
            left_half = truth_seq[:2]
            right_half = truth_seq[-2:]

        truth_seq = truth_seq.replace('I', 'B').replace('L', 'B').replace('-', '').replace('(', '').replace(')', '')

        b_seqs = [x.replace('I', 'B').replace('L', 'B') for x in b_seqs]
        y_seqs = [x.replace('I', 'B').replace('L', 'B') for x in y_seqs]

    has_left = any(
        [x == truth_seq[:len(x)] for x in b_seqs if len(x) > 1]
    ) or any(
        [truth_seq == x[:len(truth_seq)] for x in b_seqs if len(x) > 1]
    )

    has_right = any(
        [x == truth_seq[-len(x):] for x in y_seqs if len(x) > 1]
    ) or any(
        [truth_seq == x[-len(truth_seq):] for x in y_seqs if len(x) > 1]
    )

    if hybrid:
        if not has_left:
            left_half = left_half.replace('I', 'B').replace('L', 'B')
            has_left = any([left_half == x[:len(left_half)] for x in b_seqs])

        if not has_right:
            right_half = right_half.replace('I', 'B').replace('L', 'B')
            has_right = any([right_half == x[-len(right_half):] for x in y_seqs])

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



