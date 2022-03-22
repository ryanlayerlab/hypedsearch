from utils import hashable_boundaries, predicted_len
from objects import Database
import gen_spectra
from collections import defaultdict
from math import ceil
import os
import json

def write_matched_masses(filepath, matched_masses_b, matched_masses_y, kmer_set, debug):
    with open(os.path.join(filepath, "matched_masses_b.txt"),"w+") as b:
        [b.write(str(x) + ':' + json.dumps(matched_masses_b[x]) + '\n') for x in matched_masses_b.keys()]
    with open(os.path.join(filepath, "matched_masses_y.txt"),"w+") as b:
        [b.write(str(x) + ':' + json.dumps(matched_masses_y[x]) + '\n') for x in matched_masses_y.keys()]
    with open(os.path.join(filepath, "kmer_set.txt"), 'w+') as b:
        [b.write(str(x) + ':' + json.dumps(list(kmer_set[x])) + '\n') for x in kmer_set.keys()]

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

def modified_add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set,start_location,end_location,protein_number):
    for ion in 'by':
        for charge in [1, 2]:
            pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
            spec = pre_spec
            if isinstance(pre_spec,dict):
                spec = pre_spec.get('spectrum')
            for i, mz in enumerate(spec):
                start_position = start_location if ion == 'b' else end_location
                end_position = start_position + i if ion == 'b' else end_location - i
                kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
                r_d = db_dict_b if ion == 'b' else db_dict_y
                # r_d[mz].add(kmer_to_add)
                if ion == 'b':
                    r_d[mz].add((mz, protein_number, kmer_to_add, str(start_position) + '-' + str(end_position), ion, charge))
                else:
                    r_d[mz].add((mz, protein_number, kmer_to_add, str(end_position) + '-' + str(start_position), ion, charge))
                kmer_set[kmer_to_add].add(prot_name)

def make_database_set_for_protein(i,plen,max_len,prot_entry,prot_name,db_dict_b,db_dict_y,kmer_set):
    print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
    start = 1
    stop = max_len
    for j in range(start, stop):
        kmer = prot_entry.sequence[:j]
        start_position = 1
        end_position = j
        modified_add_all(kmer, prot_name, db_dict_b,db_dict_y,kmer_set, start_position, end_position, i)
    start = 0
    stop = len(prot_entry.sequence) - max_len
    for j in range(start, stop):
        kmer = prot_entry.sequence[j:j+max_len]
        start_position = j + 1
        end_position = j + max_len
        modified_add_all(kmer, prot_name, db_dict_b,db_dict_y,kmer_set,start_position, end_position, i)
    start = len(prot_entry.sequence) - max_len
    stop = len(prot_entry.sequence)
    for j in range(start, stop):
        kmer = prot_entry.sequence[j:]
        start_position = j+1
        end_position = len(prot_entry.sequence)
        modified_add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set,start_position, end_position, i)

def make_database_set_for_proteins(proteins,max_len,db_dict_b,db_dict_y,kmer_set):
    plen = len(proteins)
    for i, (prot_name, prot_entry) in enumerate(proteins):
        make_database_set_for_protein(i,plen,max_len,prot_entry,prot_name,db_dict_b,db_dict_y,kmer_set)

def modified_make_database_set(proteins: list, max_len: int):
    db_dict_b = defaultdict(set)
    db_dict_y = defaultdict(set)
    kmer_set = defaultdict(set)
    make_database_set_for_proteins(proteins,max_len,db_dict_b,db_dict_y,kmer_set)
    print('\nSorting the set of protein masses...')
    kmer_list = []
    handle_sorting_keys(db_dict_b, db_dict_y, kmer_list)
    kmer_list = sorted(kmer_list, key=lambda x: x[0])
    print('Sorting the set of protein masses done')
    return kmer_list, kmer_set

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

def modified_match_masses_per_protein(kv_prots,max_len,boundaries,kmer_set):
    extended_kv_prots = [(k, entry) for (k, v) in kv_prots for entry in v]
    kmers, kmer_set = modified_make_database_set(extended_kv_prots, max_len)
    # check_for_y_kmers(kmers)
    print("Performing Merge")
    matched_masses_b, matched_masses_y = modified_merge(kmers, boundaries)
    print("Done")
    return matched_masses_b, matched_masses_y, kmer_set

def modified_match_masses(boundaries: dict, db: Database, max_pep_len: int, debug: bool, write_path):
    max_boundary = max(boundaries.keys())
    estimated_max_len = ceil(boundaries[max_boundary][1] / 57.021464)
    max_len = min(estimated_max_len, max_pep_len)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    matched_masses_b, matched_masses_y, kmer_set = modified_match_masses_per_protein(kv_prots,max_len,boundaries,db)
    if debug:
        write_matched_masses(write_path, matched_masses_b, matched_masses_y, kmer_set, debug)
    return (matched_masses_b, matched_masses_y, kmer_set)
def reformat_kmers(kstr):
    new_list = []
    kstr = kstr.replace("[", "")
    A = kstr.rstrip().split('],')
    [new_list.append(x) for x in A]
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
def get_from_file(mb_loc, my_loc, kmer_set_loc, no_k):
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
    if no_k == False:
        with open(kmer_set_loc, 'r') as m:
            for line in m:
                line = line.replace("{", "")
                line = line.replace("}", "")
                A = line.rstrip().split(':')
                kmer_set[A[0]] = reformat_kmers(A[1])
    
    return matched_masses_b, matched_masses_y, kmer_set
