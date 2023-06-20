from sqlite import database_file
from utils import ppm_to_da
from objects import Database
import sys
import shutil
import gen_spectra
from collections import defaultdict
import os
import json
import time
from constants import ENCODED

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

# def encode_kmer(kmer):
#     # code = ''
#     # for char in kmer:
#     #     encoded = ENCODED[char]
#     #     code = code + encoded
#     # code = int(code)
#     mBytes = kmer.encode("utf-8")
#     code = int.from_bytes(mBytes, byteorder="big")
#     return code

def get_data(kmer, start, end, protein_num, ion):
    data_list = []
    for charge in [1,2]:
        # if kmer == "DPQVAQLELGG":
        #     print("here")
        mass = gen_spectra.max_mass(kmer, ion=ion, charge=charge) #1108.563316435
        ion_int = 0 if ion == 'b' else 1
        # code = encode_kmer(kmer)
        input_tuple = (mass, start, end, ion_int, charge, protein_num)
        data_list.append(input_tuple)

    return data_list
            
# def db_make_set_for_protein(i,prot,max_len, dbf, data, digest):
#     seq_len = len(prot)
#     count_max = 1000000
#     for size in range(2, max_len + 1):
#         # size -> [2, max_len]
#         for start in range(0, seq_len - size + 1):
#             end = start + size
#             kmer = prot[start:end]
#             bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
#             if not any (x in bad_chars for x in kmer):
                
#             # last_index = seq - size 6, end = start + size - 1 = 7
#             # [data.append(x) for x in get_data(kmer, start, end)]
#                 for ion in 'by':
#                     data_list = get_data(kmer, start, end, i, ion)
#                     data.extend(data_list)
#                 # insertion code
#                 if len(data) > count_max:
#                     dbf.insert(data)
#                     data.clear()
            
#     return

def db_make_set_for_protein_digest(i,prot,max_len, dbf, data, digest):
    seq_len = len(prot)
    count_max = 1000000
    for size in range(2, max_len + 1):
        # size -> [2, max_len]
        for start in range(0, seq_len - size + 1):
            end = start + size
            kmer = prot[start:end]
            if kmer[0] in digest[0] or digest[0] == ['-'] or (start > 0 and prot[start-1] in digest[1]): #cuts to the left
                bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'b')
                    data.extend(data_list)
                    # insertion code
                    if len(data) > count_max:
                        dbf.insert(data)
                        data.clear()
            if kmer[-1] in digest[1] or digest[1] == ['-'] or (end < seq_len and prot[end] in digest[0]): #cuts to the right
                bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'y')
                    data.extend(data_list)
                    # insertion code
                    if len(data) > count_max:
                        dbf.insert(data)
                        data.clear()
            
    return

def db_make_database_set_for_proteins(proteins,max_len,dbf,digest):
    plen = len(proteins)
    last_percent = 0
    data = []
    
    for i, (_, prot_entry) in enumerate(proteins):
        percent = int((i+1) * 100 / plen)
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        if percent != last_percent:
            # print(f'\rInserting {percent}%', end='')
            last_percent = percent
            free = shutil.disk_usage('/')[2]
            free = free/(1024**3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB" )
                sys.exit(1)
        db_make_set_for_protein_digest(i,prot_entry,max_len, dbf, data, digest)
        
    if len(data) != 0:
        dbf.insert(data)
        
def modified_make_database_set(proteins: list, max_len: int, dbf, digest):
    
    print("\nBeginning Insertions")
    start = time.time()
    db_make_database_set_for_proteins(proteins,max_len,dbf,digest)
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

def modified_match_masses(input_masses: list, db: Database, max_len: int, ppm_tolerance, b_prec, y_prec):
    dbf = database_file(max_len, False)
    matched_masses_b, matched_masses_y = dict(), dict()
    
    start = time.time()
    for input_mass in input_masses:
        tol = ppm_to_da(input_mass, ppm_tolerance)
        matched_masses_b[input_mass], matched_masses_y[input_mass] = dbf.query_mass(input_mass, tol)
        
    tol = ppm_to_da(b_prec, ppm_tolerance)
    matched_masses_b[b_prec], _ = dbf.query_mass(b_prec, tol)
    tol = ppm_to_da(y_prec, ppm_tolerance)
    _, matched_masses_y[y_prec] = dbf.query_mass(y_prec, tol)
        
    end = time.time() - start
    with open('Timing_data.txt', 'w') as t:
        t.write("Queries took:" + '\t' + str(end) + "\n")       

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
