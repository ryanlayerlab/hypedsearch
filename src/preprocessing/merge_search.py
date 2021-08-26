from utils import hashable_boundaries, predicted_len
from objects import Database
import gen_spectra
from collections import defaultdict
from typing import Iterable
from math import ceil
import array as arr

BATCH_SIZE = 300

def merge(mz_s: Iterable, indices: Iterable, kmers: Iterable, boundaries: Iterable):
    boundry_index, mass_index = 0, 0
    matched_masses = defaultdict(list)
    while boundry_index < len(boundaries) and mass_index < len(mz_s):
       if boundaries[boundry_index][0] <= mz_s[mass_index] <= boundaries[boundry_index][1]:
           matched_masses[hashable_boundaries(boundaries[boundry_index])] += kmers[indices[mass_index - 1]:indices[mass_index]]
           mass_index += 1
       elif mz_s[mass_index] > boundaries[boundry_index][1]:
           boundry_index += 1
       elif mz_s[mass_index] < boundaries[boundry_index][0]:
           mass_index += 1
    return matched_masses

def add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set):
    for ion in 'by':
        for charge in [1, 2]:
            pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
            spec = pre_spec
            if isinstance(pre_spec,dict):
                spec = pre_spec.get('spectrum')
            for i, mz in enumerate(spec):
                kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
                r_d = db_dict_b if ion == 'b' else db_dict_y
                r_d[mz].add(kmer_to_add)
                kmer_set[kmer_to_add].append(prot_name)

def make_database_set_for_protein(i,plen,max_len,prot_entry,prot_name,db_dict_b,db_dict_y,kmer_set):
    print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
    start = 1
    stop = max_len
    for j in range(start, stop):
        kmer = prot_entry.sequence[:j]
        add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set)
    start = 0
    stop = len(prot_entry.sequence) - max_len
    for j in range(start,stop):
        kmer = prot_entry.sequence[j:j+max_len]
        add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set)
    start = len(prot_entry.sequence) - max_len
    stop = len(prot_entry.sequence)
    for j in range(start, stop):
        kmer = prot_entry.sequence[j:]
        add_all(kmer, prot_name,db_dict_b,db_dict_y,kmer_set)

def sort_masses_in_sorted_keys_b(db_dict_b,mz,db_list_b,index_list_b,kmer_list_b):
    kmers = db_dict_b[mz]
    db_list_b.append(mz)
    offset = 0 if not len(index_list_b) else index_list_b[-1]
    index_list_b.append(len(kmers) + offset)
    kmer_list_b += kmers

def sort_masses_in_sorted_keys_y(db_dict_y,mz,db_list_y,index_list_y,kmer_list_y):
    kmers = db_dict_y[mz]
    db_list_y.append(mz)
    offset = 0 if not len(index_list_y) else index_list_y[-1]
    index_list_y.append(len(kmers) + offset)
    kmer_list_y += kmers

def handle_sorting_keys_b(db_dict_b,db_list_b,index_list_b,kmer_list_b):
    sorted_keys = sorted(db_dict_b.keys())
    for mz in sorted_keys:
        sort_masses_in_sorted_keys_b(db_dict_b,mz,db_list_b,index_list_b,kmer_list_b)

def handle_sorting_keys_y(db_dict_y,db_list_y,index_list_y,kmer_list_y):
    sorted_keys = sorted(db_dict_y.keys())
    for mz in sorted_keys:
        sort_masses_in_sorted_keys_y(db_dict_y,mz,db_list_y,index_list_y,kmer_list_y)

def make_database_set_for_proteins(proteins,max_len,db_dict_b,db_dict_y,kmer_set):
    plen = len(proteins)
    for i, (prot_name, prot_entry) in enumerate(proteins):
        make_database_set_for_protein(i,plen,max_len,prot_entry,prot_name,db_dict_b,db_dict_y,kmer_set)

def make_database_set(proteins: list, max_len: int):
    db_dict_b = defaultdict(set)
    db_dict_y = defaultdict(set)
    kmer_set = defaultdict(list)
    make_database_set_for_proteins(proteins,max_len,db_dict_b,db_dict_y,kmer_set)
    print('\nSorting the set of protein masses...')
    db_list_b, index_list_b, kmer_list_b = arr.array('f'), arr.array('i'), []
    db_list_y, index_list_y, kmer_list_y = arr.array('f'), arr.array('i'), []
    handle_sorting_keys_b(db_dict_b,db_list_b,index_list_b,kmer_list_b)
    handle_sorting_keys_y(db_dict_y,db_list_y,index_list_y,kmer_list_y)
    print('Sorting the set of protein masses done')
    return db_list_b, index_list_b, kmer_list_b, db_list_y, index_list_y, kmer_list_y, kmer_set

def add_matched_to_matched_set(matched_masses_b_batch,matched_masses_b,kmer_set,batch_kmer_set,matched_masses_y_batch,matched_masses_y):
    for k, v in matched_masses_b_batch.items():
        matched_masses_b[k] += v 
        for kmer in v:
            kmer_set[kmer] += batch_kmer_set[kmer]
    for k, v in matched_masses_y_batch.items():
        matched_masses_y[k] += v 
        for kmer in v:
            kmer_set[kmer] += batch_kmer_set[kmer]

def match_masses_per_protein(batch_num,num_batches,max_len,batch_set,spectra_boundaries,kmer_set, matched_masses_b,matched_masses_y):
    print(f'On batch {batch_num + 1}/{num_batches}\n', end='')
    extended_batch_set = [(k, entry) for (k, v) in batch_set for entry in v]
    batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = make_database_set(extended_batch_set, max_len)
    matched_masses_b_batch = merge(batch_b_list, index_list_b, batch_kmer_b, spectra_boundaries)
    matched_masses_y_batch = merge(batch_y_list, index_list_y, batch_kmer_y, spectra_boundaries)
    add_matched_to_matched_set(matched_masses_b_batch,matched_masses_b,kmer_set,batch_kmer_set,matched_masses_y_batch,matched_masses_y)

def match_masses_for_all_proteins(batched_prots,num_batches,max_len,spectra_boundaries,kmer_set,matched_masses_b,matched_masses_y):
    for batch_num, batch_set in enumerate(batched_prots):
        match_masses_per_protein(batch_num,num_batches,max_len,batch_set,spectra_boundaries,kmer_set, matched_masses_b,matched_masses_y)

def match_masses(spectra_boundaries: list, db: Database, max_pep_len: int = 30):
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)
    estimated_max_len = ceil(spectra_boundaries[-1][1] / 57.021464)
    max_len = min(estimated_max_len, max_pep_len)
    num_batches = ceil(len(db.proteins) / BATCH_SIZE)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]
    match_masses_for_all_proteins(batched_prots,num_batches,max_len,spectra_boundaries,kmer_set,matched_masses_b,matched_masses_y)
    db = db._replace(kmers=kmer_set)
    return (matched_masses_b, matched_masses_y, db)