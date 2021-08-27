import os, gzip, shutil, copy, math
from typing import Iterable, Any
from itertools import product
import numpy as np
from collections import namedtuple

from objects import Spectrum
import gen_spectra

import math
import re

HYBRID_ALIGNMENT_PATTERN = re.compile(r'[-\(\)]')

def file_exists(file_name: str) -> bool:
    '''Determine if a file exists

    :param file_name: Path to the file in question
    :type file_name: str
    
    :returns: True if the file exists
    :rtype: bool
    '''
    
    return os.path.isfile(file_name)

def make_valid_dir_string(dir_path: str) -> str:
    '''
    Add os separator character to end of directory string to make valid directory path

    :param dir_path: Name of directory to check 
    :type dir_path: str

    :returns: Corrected directory path
    :rtype: str
    '''

    return dir_path + os.path.sep if os.path.sep != dir_path[-1] else dir_path

def make_dir(dir_path: str) -> bool:
    '''Check directory path for existing directory or make one with the name given. 
    NOTE: this is not recursive, only 1 level of directory will be created

    :param dir_path: Full path of directory to create
    :type dir_path: str

    :returns: True if successful
    :rtype: bool
    '''

    dir_path = make_valid_dir_string(dir_path)
    try:
        if not os.path.exists(dir_path): 
            os.makedirs(dir_path)
        return True 
    except:
        return False

def make_valid_text_file(file_name: str) -> str:
    '''Ensure some string path has .txt appended to it for appropriate .txt extension

    :param file_name: File name to validate for txt file type
    :type file_name: str

    :returns: Name with the .txt extension
    :rtype: str
    '''

    file_name = file_name + '.txt' if '.txt' not in file_name else file_name
    return file_name

def make_valid_json_file(file_name: str) -> str:
    '''Ensure some string path has .json appended to it for appropriate .json extension

    :param file_name: File name to validate for json file type
    :type file_name: str

    :returns: Name with the .json extension
    :rtype: str
    '''

    file_name = file_name + '.json' if '.json' not in file_name else file_name
    return file_name

def make_valid_csv_file(file_name: str) -> str:
    '''Ensure some string path has .csv appended to it for appropriate .csv extension

    :param file_name: File name to validate for csv file type
    :type file_name: str

    :returns: Name with the .csv extension
    :rtype: str
    '''

    file_name = file_name + '.csv' if '.csv' not in file_name else file_name
    return file_name

def make_valid_fasta_file(file_name: str) -> str:
    '''Ensure some string path has .fasta appended to it for appropriate .fasta extension

    :param file_name: File name to validate for fasta file type
    :type file_name: str

    :returns: Name with the .fasta extension
    :rtype: str
    '''

    file_name = file_name + '.fasta' if '.fasta' not in file_name else file_name
    return file_name

def is_json(file: str) -> bool:
    '''Determine if a file is a json file

    :param file: File name of file in question
    :type file: str

    :returns: True if is a json file
    :rtype: bool
    '''

    return True if '.json' in file else False

def is_fasta(file: str) -> bool:
    '''Determine if a file is a fasta file 

    :param file: File name of the file in question
    :type file: str

    :returns: True if it is a fasta file
    :rtype: str 
    '''

    return True if '.fasta' in file else False

def is_dir(dir_path: str) -> bool:
    '''Determine if a path is a valid path to a directory

    :param dir_path: Full path tl the director
    :type dir_path: str

    :returns: True if the directory exists
    :rtype: bool
    '''

    return os.path.isdir(dir_path)

def is_file(file: str) -> bool:
    '''Determine if a file exists

    :param file: Full path to the file
    :type file: str
    
    :returns: True if the file exists
    :rtype: bool
    '''

    return os.path.isfile(file)

def all_perms_of_s(s: str, keyletters: str) -> list:
    '''Find all permutations of a string that has values 'keyletters' in them

    :param s: The string to permutate
    :type s: str
    :param keyletters: The characters to make permutations with
    :type keyletters: str

    :returns: All permutations of s
    :rtype: list

    :Example:

    >>> all_perms_of_s('LMNOP', 'LI')
    >>> ['LMNOP', 'IMNOP']
    '''

    # Convert input string into a list so we can easily substitute letters
    seq = list(s)
    
    perms = []

    # Find indices of key letters in seq
    indices = [ i for i, c in enumerate(seq) if c in keyletters ]

    # Generate key letter combinations & place them into the list
    for t in product(keyletters, repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        perms.append(''.join(seq))
    return perms

def ppm_to_da(mass: float, ppm_tolerance: float) -> float:
    return abs((ppm_tolerance / 1000000)*mass)

def make_sparse_array(spectrum: list, width: float, value=50) -> np.ndarray:
    '''Make a spectrum (a list of floats) into a sparsely populated array for xcorr 
    calculation. Indices are calculated by
    
    idx = int(m/w), m is mass, w is bin width
    
    width is the tolerance in Da to allow when calculating scores. All peaks with some value
    are given a new value of 50.
    
    :param spectrum: Floating point mass values of peaks
    :type spectrum: list
    :param width: Mass tolerance for bin width
    :type width: foat
    :param value: Value to put in a bin where a mass is found
        (default is 50)
    :type value: number

    :returns: Sparesly populated value-hot array
    :rtype: numpy.ndarray
    '''

    # find the largest mass and make that the length of the array
    list_size = int(max(spectrum)//width) + 1
    
    sparse = np.zeros(list_size)
    
    # populate sparse at the index for each mass
    for m in spectrum:
        sparse[int(m // width)] = value

    return sparse

def overlap_intervals(intervals: list) -> list:
    '''Take a list of intervals and turn it into a smaller list by finding any 
    overlapping intervals and making it a larger interval

    :param intervals: Intervals (in the form of lists [lower_bound, upper_bound]). Both ends are inclusive
    :type intervals: list

    :returns: Overlapped intervals of [lower_bound, upper_bound]
    :rtype: list 
    '''

    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged  

def to_percent(index, total):
    return int(100 * (index)/total)

def predicted_len(precursor_mass: float, precursor_charge: int) -> int:
    '''The predicted length of a spectrum based on its maximum mass

    :param precursor_mass: The maximum mass of the sequence
    :type precursor_mass: float
    :param precursor_charge: The charge of the observed precusor mass
    :type precursor_charge: int

    :returns: Predicted sequence length
    :rtype: int
    '''

    return math.ceil(precursor_mass / gen_spectra.get_precursor('G', precursor_charge))

def predicted_len_precursor(spectrum: Spectrum, sequence: str) -> int:
    '''Make a prediction of the peptide length give a spectrum and the current 
    sequence. 

    :param spectrum: The observed spectrum 
    :type spectrum: Spectrum
    :param sequence: The current alignment made
    :type sequence: str

    :returns: Predicted length of a full alignment
    :rtype: int
    '''

    # first get the theoretical precursor of the sequence
    theoretical_prec = gen_spectra.get_precursor(sequence, spectrum.precursor_charge)

    # now we can estimate length   spec/seq = real/theory
    estimated_len = math.ceil(len(sequence) * (spectrum.precursor_mass / theoretical_prec))

    return estimated_len

def hashable_boundaries(boundaries: list) -> str:
    '''Turn a lower and upper bound into a string in order 
    to hash

    :param boundaries: A list of lists where each internal list is [lower_bound, upper_bound]
    :type boundaries: list

    :returns: A string of the lower and upper bounds connected that looks like <lower_bound>-<upper_bound>
    :rtype: str
    '''
    if (len(boundaries) == 2):
        return '-'.join([str(x) for x in boundaries])
    else:
        return None

def cosine_similarity(a: list, b: list) -> float:
    '''Calculate the cosine similarity of two vectors
    
    :param a: First vector 
    :type a: list
    :param b: Second vector 
    :type b: list

    :returns: The cosine similarity of the two vectors
    :rtype: float
    '''

    # pad with zeros if not the same size
    if len(a) > len(b):
        b = list(b) + [0 for _ in range(len(a) - len(b))]
        
    if len(b) > len(a):
        a = list(a) + [0 for _ in range(len(b) - len(a))]

    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def split_hybrid(sequence: str) -> (str, str):
    '''Split a hybrid sequence into it's left and right components
    
    :param sequence: hybrid sequence with special characters [() -]
    :type sequence: str
    
    :returns: left subsequence, right subsequence
    :rtype: (str, str)
    '''
    if '-' in sequence:
        return (sequence.split('-')[0], sequence.split('-')[1])
    
    else:
        left = sequence.split(')')[0].replace('(', '')
        right = sequence.split('(')[1].replace(')', '')
        return (left, right)

def DEV_contains_truth_parts(truth_seq: str, hybrid: bool, b_seqs: list, y_seqs: list) -> bool:
    '''DEVELOPMENT FUNCTION ONLY

    Determines if a set of b and y sequences can potentially create the truth. If so, 
    True is returned, otherwise False

    :param truth_seq: The "truth" sequence or the sequence we want to try and find for this spectrum
    :type truth_seq: str
    :param hybrid: Is the alignment supposed to be a hybrid sequence
    :type hybrid: bool 
    :param b_seqs: k-mers identified from the b-ion score
    :type b_seqs: list
    :param y_seqs: k-mers identified from the y-ion score
    :type y_seqs: list

    :returns: Whether or not the "true" or desired sequence could be found from the b and y sequences
    :rtype: bool
    '''
    
    has_left = False
    has_right = False

    left_half = ''
    right_half = ''

    # if the sequence is hybrid, replace all I and L with B
    if hybrid:
        
        # split into left and right halves for later 
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

    # see if any of the b seqs match the first bit of the sequence
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

    # if its a hybrid, both left and right must be true, otherwise just one will do 
    if hybrid:
        
        # if we dont' have complete matches, thats ok. If we can get the halves 
        # of the "truth" in some longer b or y seq, thats good enough
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
    '''DEVELOPMENT FUNCTION ONLY

    Determines of the truth sequence is held in the list of sequences. If so, 
    True is returned, otherwise False

    :param truth_seq: The "truth" sequence or the sequence we want to try and find for this spectrum
    :type truth_seq: str
    :param hybrid: Is the alignment supposed to be a hybrid sequence
    :type hybrid: bool
    :param seqs: The sequences to look through
    :type seqs: list

    :returns: True if the "truth" sequence is found in the list
    :rtype: bool
    '''

    # if hybrid, replace all L and I with B, and also remove any - or ()
    if hybrid:
        truth_seq = truth_seq \
                        .replace('I', 'B') \
                        .replace('L', 'B') \
                        .replace('-', '') \
                        .replace('(', '') \
                        .replace(')', '')

        seqs = [
            x \
                .replace('I', 'B') \
                .replace('L', 'B') \
                .replace('-', '') \
                .replace('(', '') \
                .replace(')', '')
            for x in seqs
        ]

    contains_exact = any([x == truth_seq for x in seqs])

    return contains_exact