from xmlrpc.client import Boolean

from pytz import country_names
from sqlite import database_file
from utils import hashable_boundaries, ppm_to_da, predicted_len
from objects import Database
import sys
import shutil
import gen_spectra
from collections import defaultdict
from math import ceil
import os
import json
import time

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

# def modified_add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set,start_location,end_location,protein_number):
#     for ion in 'by':
#         for charge in [1, 2]:
#             pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
#             spec = pre_spec
#             if isinstance(pre_spec,dict):
#                 spec = pre_spec.get('spectrum')
#             for i, mz in enumerate(spec):
#                 start_position = start_location if ion == 'b' else end_location
#                 end_position = start_position + i if ion == 'b' else end_location - i
#                 kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
#                 r_d = db_dict_b if ion == 'b' else db_dict_y
#                 # r_d[mz].add(kmer_to_add)
#                 if ion == 'b':
#                     r_d[mz].add((mz, protein_number, kmer_to_add, str(start_position) + '-' + str(end_position), ion, charge))
#                 else:
#                     r_d[mz].add((mz, protein_number, kmer_to_add, str(end_position) + '-' + str(start_position), ion, charge))
#                 kmer_set[kmer_to_add].add(prot_name)

# def make_database_set_for_protein(i,plen,max_len,prot_entry,prot_name,db_dict_b,db_dict_y,kmer_set):
#     print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
#     start = 0
#     stop = max_len
#     for j in range(start, stop):
#         kmer = prot_entry.sequence[:j+1]
#         start_position = start
#         end_position = j
#         modified_add_all(kmer, prot_name, db_dict_b,db_dict_y,kmer_set, start_position, end_position, i)
#     start = 0
#     stop = len(prot_entry.sequence) - max_len
#     for j in range(start, stop):
#         kmer = prot_entry.sequence[j:j+max_len]
#         start_position = j
#         end_position = j + max_len - 1
#         modified_add_all(kmer, prot_name, db_dict_b,db_dict_y,kmer_set,start_position, end_position, i)
#     start = len(prot_entry.sequence) - max_len
#     stop = len(prot_entry.sequence)
#     for j in range(start, stop):
#         kmer = prot_entry.sequence[j:]
#         start_position = j
#         end_position = len(prot_entry.sequence)
#         modified_add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set,start_position, end_position, i)

def get_data(kmer, start, end, protein_num):
    data_list = []
    for ion in 'by':
        for charge in [1,2]:
            mass = gen_spectra.max_mass(kmer, ion=ion, charge=charge)
            ion_int = 0 if ion == 'b' else 1
            input_tuple = (mass, start, end, ion_int, charge, protein_num)
            data_list.append(input_tuple)

    return data_list
            
def db_make_set_for_protein(i,prot,max_len, dbf, data):
    seq_len = len(prot)
    count_max = 1000000
    for size in range(2, max_len + 1):
        # size -> [2, max_len]
        for start in range(0, seq_len - size + 1):
            end = start + size
            kmer = prot[start:end]
            # last_index = seq - size 6, end = start + size - 1 = 7
            # [data.append(x) for x in get_data(kmer, start, end)]
            data_list = get_data(kmer, start, end, i)
            data.extend(data_list)
            # insertion code
            if len(data) > count_max:
                dbf.insert(data)
                data.clear()
            
    return

def db_make_database_set_for_proteins(proteins,max_len,dbf):
    plen = len(proteins)
    last_percent = 0
    data = []
    
    for i, (_, prot_entry) in enumerate(proteins):
        percent = int((i+1) * 100 / plen)
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        if i == 256:
            continue
        if percent != last_percent:
            # print(f'\rInserting {percent}%', end='')
            last_percent = percent
            free = shutil.disk_usage('/')[2]
            free = free/(1024**3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB" )
                sys.exit(1)
        db_make_set_for_protein(i,prot_entry,max_len, dbf, data)
        
    if len(data) != 0:
        dbf.insert(data)
        
def modified_make_database_set(proteins: list, max_len: int, dbf):    
    # final data structure: dict-> called mass_kmers
    # input_masses -> list of tuples
    # mass_kmers[mass] -> list of tuples
    # db.query_mass(mass) -> list of tuples
    # for input_mass in input_masses:
    #   kmers = db.query_mass(input_mass) 
    #   if len(kmers) != 0:
    #       we have the list
    
    
    #kmer -> list of proteins where it occurs
    #kmer list: list of tuples (kmer, mass, location, ion, charge)
    print("\nBeginning Insertions")
    start = time.time()
    db_make_database_set_for_proteins(proteins,max_len, dbf)
    duration = time.time() - start
    print("Insertion took: ", duration)
    # db.read() #Only for debugging
    print('\nIndexing the set of kmers based on mass, ion')
    dbf.index_ion_mass()
    print('\nIndexing the set of kmers based on protein, start position, end position')
    dbf.index_ion_mass_b()
    print('\nIndexing the set of kmers based on protein, end position, start position')
    dbf.index_ion_mass_y()
    print('Done making database')
    
    
    # print('\nSorting the set of protein masses...')
    # kmer_list = []
    # # handle_sorting_keys(db_dict_b, db_dict_y, kmer_list)
    # kmer_list = sorted(kmer_list, key=lambda x: x[0])
    # print('Sorting the set of protein masses done')
    return

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

def modified_match_masses(input_masses: list, db: Database, max_len: int, ppm_tolerance, make_new):
    # max_boundary = max(boundaries.keys())
    # estimated_max_len = ceil(boundaries[max_boundary][1] / 57.021464)
    # max_len = min(estimated_max_len, max_pep_len)
    
    dbf = database_file(max_len, make_new)
    if make_new:
        kv_prots = [(k, v) for k, v in db.proteins]    
        # extended_kv_prots = [(k, entry) for (k, v) in kv_prots for entry in v]
        modified_make_database_set(kv_prots, max_len, dbf)
    
    matched_masses_b, matched_masses_y = dict(), dict()
    
    start = time.time()
    for input_mass in input_masses:
        tol = ppm_to_da(input_mass, ppm_tolerance)
        # ion_int = 0 if ion == 'b' else 1
        matched_masses_b[input_mass], matched_masses_y[input_mass] = dbf.query_mass(input_mass, tol) #same place: location start, protein_num
        # print(input_mass, 'b', matched_masses_b[input_mass])
        # print(input_mass, 'y', matched_masses_y[input_mass])
        
    end = time.time() - start
    with open('Timing_data.txt', 'w') as t:
        t.write("Queries took:" + '\t' + str(end) + "\n")       
        # if debug:
        # write_matched_masses(write_path, matched_masses_b, matched_masses_y, kmer_set, debug)

    return matched_masses_b, matched_masses_y
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
