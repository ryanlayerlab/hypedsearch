from lookups.sqlite_database import Sqllite_Database
from lookups.utils import ppm_to_da
from lookups.objects import FastaDatabase
import sys
import shutil
import computational_pipeline.gen_spectra
from collections import defaultdict
import os
import json
import time
from lookups.constants import ENCODED

def write_matched_masses(filepath, matched_masses_b, matched_masses_y, kmer_set, debug):
    with open(os.path.join(filepath, "matched_masses_b.txt"),"w+") as b:
        [b.write(str(x) + ':' + json.dumps(matched_masses_b[x]) + '\n') for x in matched_masses_b.keys()]
    with open(os.path.join(filepath, "matched_masses_y.txt"),"w+") as b:
        [b.write(str(x) + ':' + json.dumps(matched_masses_y[x]) + '\n') for x in matched_masses_y.keys()]
    with open(os.path.join(filepath, "kmer_set.txt"), 'w+') as b:
        [b.write(x + ':' + json.dumps(list(kmer_set[x])) + '\n') for x in kmer_set.keys()]

def modified_sort_masses_in_sorted_keys_b(db_dict_b,mz,kmer_list_b):
    kmers = db_dict_b[mz]
    kmer_list_b += kmers

def modified_sort_masses_in_sorted_keys_y(db_dict_y,mz,kmer_list_y):
    kmers = db_dict_y[mz]
    kmer_list_y += kmers

def handle_sorting_keys(db_dict_b, db_dict_y, kmer_list):
    sorted_b_keys = sorted(db_dict_b.keys())
    sorted_y_keys = sorted(db_dict_y.keys())
    for mz in sorted_b_keys:
        modified_sort_masses_in_sorted_keys_b(db_dict_b,mz,kmer_list)
    for mz in sorted_y_keys:
        modified_sort_masses_in_sorted_keys_y(db_dict_y,mz,kmer_list)

def in_bounds(int1, interval):
    if int1 >= interval[0] and int1 <= interval[1]:
        return True
    else:
        return False

def modified_merge(kmers, boundaries: dict):
    matched_masses_b, matched_masses_y = defaultdict(list), defaultdict(list)
    mz_mapping = dict()
    for i,mz in enumerate(boundaries):
        mz_mapping[i] = []
        mz_mapping[i].append(mz)
        [mz_mapping[i].append(x) for x in boundaries[mz]]
    boundary_index, kmer_index, starting_point = 0,0,0
    while (boundary_index < len(boundaries)) and (kmer_index < len(kmers)):
        target_kmer = kmers[kmer_index]
        target_key = mz_mapping[boundary_index][0]
        target_boundary = mz_mapping[boundary_index][1:]
        if in_bounds(target_kmer[0], target_boundary):
            if target_kmer[4] == 'b':
                matched_masses_b[target_key].append(target_kmer)
                kmer_index = kmer_index + 1

            if target_kmer[4] == 'y':
                matched_masses_y[target_key].append(target_kmer)
                kmer_index = kmer_index + 1
        elif target_kmer[0] < target_boundary[0]:
            kmer_index = kmer_index + 1
            starting_point = starting_point + 1
        else:                                           
            boundary_index = boundary_index + 1
            kmer_index = starting_point
    return matched_masses_b, matched_masses_y

def get_all_matched_rows(input_masses, sqllite_database, ppm_tolerance):
    all_b_rows = []
    all_y_rows = []
    for input_mass in input_masses:
        tolerance = ppm_to_da(input_mass, ppm_tolerance)
        b_rows, y_rows = sqllite_database.query_mass_kmers(input_mass, tolerance)
        all_b_rows.extend(b_rows)
        all_y_rows.extend(y_rows)
    return all_b_rows, all_y_rows

def reformat_kmers(kstr):
    new_list = []
    kstr = kstr.replace("[", "")
    kstr = kstr.replace("]", "")
    kstr = kstr.replace(" ", "")
    kstr = kstr.replace("\"", "")
    A = kstr.rstrip().split(',')
    [new_list.append(str(x)) for x in A]
    return new_list
def reformat_hits(cstr):
    new_list = []
    cstr = cstr.replace("[", "")
    cstr = cstr.replace(" ", "")
    cstr = cstr[:-2]
    A = cstr.rstrip().split('],')
    for block in A:
        sublist = []
        B = block.rstrip().split(',')
        sublist.append(float(B[0]))
        sublist.append(int(B[1]))
        sublist.append(B[2][1:-1])
        sublist.append(B[3][1:-1])
        sublist.append(B[4][1:-1])
        sublist.append(int(B[5]))
        new_list.append(sublist)
    return new_list
def get_from_file(mb_loc, my_loc, kmer_set_loc, kmer):
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(), defaultdict(), defaultdict()
    with open(mb_loc, 'r') as m:
        for line in m:
            line = line.replace("{", "")
            line = line.replace("}", "")
            A = line.rstrip().split(':')
            matched_masses_b[float(A[0])] = reformat_hits(A[1])
    with open(my_loc, 'r') as m:
        for line in m:
            line = line.replace("{", "")
            line = line.replace("}", "")
            A = line.rstrip().split(':')
            matched_masses_y[float(A[0])] = reformat_hits(A[1])      
    if kmer == True:
        with open(kmer_set_loc, 'r') as m:
            for line in m:
                line = line.replace("{", "")
                line = line.replace("}", "")
                A = line.rstrip().split(':')
                kmer_set[A[0]] = reformat_kmers(A[1])
    
    return matched_masses_b, matched_masses_y, kmer_set
