from computational_pipeline.sqlite import database_file
from lookups.utils import ppm_to_da
from lookups.objects import Database
import computational_pipeline.gen_spectra
from lookups.constants import ENCODED
from collections import defaultdict
import sys, shutil, os, json

def write_matched_masses(filepath, matched_masses_b, matched_masses_y, kmer_set, debug):
    def write_to_file(filename, dictionary):
        with open(os.path.join(filepath, filename), 'w+') as file:
            for key, val in dictionary.items():
                file.write(f'{key}:{json.dumps(val)}\n')
    write_to_file('matched_masses_b.txt', matched_masses_b)
    write_to_file('matched_masses_y.txt', matched_masses_y)
    write_to_file('kmer_set.txt', kmer_set)

def modified_sort_masses_in_sorted_keys_b(db_dict_b,mz,kmer_list_b):
    kmers = db_dict_b[mz]
    kmer_list_b += kmers

def modified_sort_masses_in_sorted_keys_y(db_dict_y,mz,kmer_list_y):
    kmers = db_dict_y[mz]
    kmer_list_y += kmers

def handle_sorting_keys(db_dict_b, db_dict_y, kmer_list):
    sorted_b_keys = sorted(db_dict_b)
    sorted_y_keys = sorted(db_dict_y)
    for mz in sorted_b_keys:
        modified_sort_masses_in_sorted_keys_b(db_dict_b, mz, kmer_list)
    for mz in sorted_y_keys:
        modified_sort_masses_in_sorted_keys_y(db_dict_y, mz, kmer_list)

def get_data(kmer, start, end, protein_num, ion):
    data_list = []
    ion_int = int(ion != 'b')
    for charge in [1,2]:
        mass = computational_pipeline.gen_spectra.max_mass(kmer, ion=ion, charge=charge) #1108.563316435
        input_tuple = (mass, start, end, ion_int, charge, protein_num)
        data_list.append(input_tuple)
    return data_list
            
def db_make_set_for_protein_digest(i, prot, max_len, dbf, data, digest):
    seq_len = len(prot)
    COUNT_MAX = 1000000
    for size in range(2, max_len + 1):
        # size -> [2, max_len]
        for start in range(0, seq_len - size + 1):
            end = start + size
            kmer = prot[start:end]
            bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
            if kmer[0] in digest[0] or digest[0] == ['-'] or (start > 0 and prot[start-1] in digest[1]): #cuts to the left
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'b')
                    data.extend(data_list)
                    # insertion code
                    if len(data) > COUNT_MAX:
                        dbf.insert(data)
                        data.clear()
            if kmer[-1] in digest[1] or digest[1] == ['-'] or (end < seq_len and prot[end] in digest[0]): #cuts to the right
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'y')
                    data.extend(data_list)
                    # insertion code
                    if len(data) > COUNT_MAX:
                        dbf.insert(data)
                        data.clear()

def db_make_database_set_for_proteins(proteins,max_len,dbf,digest):
    plen = len(proteins)
    prev_percent = 0
    data = []
    
    for i, (_, prot_entry) in enumerate(proteins):
        percent = int((i+1) * 100 / plen)
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        if percent != prev_percent:
            prev_percent = percent
            space_available = shutil.disk_usage('/')[2] / (1024**3)
            if space_available < 10:
                print("\nUsed too much space, Space available =", space_available, "GB" )
                sys.exit(1)
        db_make_set_for_protein_digest(i, prot_entry, max_len, dbf, data, digest)
        
    if data:
        dbf.insert(data)
        
def modified_make_database_set(proteins: list, max_len: int, dbf, digest):
    db_make_database_set_for_proteins(proteins,max_len,dbf,digest)
    dbf.index_ion_mass()
    dbf.index_ion_mass_b()
    dbf.index_ion_mass_y()

def in_bounds(int_var, interval):
    return int_var >= interval[0] and int_var <= interval[1]

def modified_merge(kmers, boundaries: dict):
    matched_masses_b, matched_masses_y = defaultdict(list), defaultdict(list)
    mz_mapping = dict()
    for i, mz in enumerate(boundaries):
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
    
    for input_mass in input_masses:
        tol = ppm_to_da(input_mass, ppm_tolerance)
        matched_masses_b[input_mass], matched_masses_y[input_mass] = dbf.query_mass(input_mass, tol)
        
    tol = ppm_to_da(b_prec, ppm_tolerance)
    matched_masses_b[b_prec], _ = dbf.query_mass(b_prec, tol)
    tol = ppm_to_da(y_prec, ppm_tolerance)
    _, matched_masses_y[y_prec] = dbf.query_mass(y_prec, tol)
        

    return matched_masses_b, matched_masses_y

def reformat_kmers(kstr):
    kstr = kstr.replace("[", "")
    kstr = kstr.replace("]", "")
    kstr = kstr.replace(" ", "")
    kstr = kstr.replace("\"", "")
    new_list = kstr.rstrip().split(',')
    return new_list

def reformat_hits(cstr):
    new_list = []
    cstr = cstr.replace("[", "")
    cstr = cstr.replace(" ", "")
    cstr = cstr[:-2]
    A = cstr.rstrip().split('],')
    for block in A:
        B = block.rstrip().split(',')
        sublist = [
            float(B[0]), int(B[1]), B[2][1:-1], 
            B[3][1:-1], B[4][1:-1], int(B[5])
        ]
        new_list.append(sublist)
    return new_list

def get_from_file(mb_loc, my_loc, kmer_set_loc, kmer):
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(), defaultdict(), defaultdict()
    def do_thing(dictionary, filename, func):
        with open(filename, 'r') as file: 
            for line in file: 
                line = line.replace('{', '').replace('}', '')
                split_line = line.rstrip().split(':')
                dictionary[float[split_line[0]]] = func(split_line[1])
                
    do_thing(matched_masses_b, mb_loc, reformat_hits)
    do_thing(matched_masses_y, my_loc, reformat_hits)
    if kmer: 
        do_thing(kmer_set, kmer_set_loc, reformat_kmers)
    return matched_masses_b, matched_masses_y, kmer_set
