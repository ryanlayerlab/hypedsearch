from utils import hashable_boundaries, predicted_len
from objects import Database
import gen_spectra
from collections import defaultdict
from typing import Iterable
from math import ceil
import array as arr
import requests
import json

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
                kmer_set[kmer_to_add].append(prot_name)

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
    kmer_set = defaultdict(list)
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
    #Goal: b and y dictionaries mapping mz values to lists of kmers that have a mass within the tolerance
    # kmers = make_database_set(db.proteins, max_len)
    mz_mapping = dict()
    for i,mz in enumerate(boundaries):
        mz_mapping[i] = boundaries[mz]
    boundary_index, kmer_index, starting_point = 0,0,0
    while (boundary_index < len(boundaries)) and (kmer_index < len(kmers)):
        #idea is to increment kmer index when mass is too small for boundaries[0] and then stop when mass is too big for boundaries[1]
        target_kmer = kmers[kmer_index]
        target_boundary = mz_mapping[boundary_index]
        if in_bounds(target_kmer[0], target_boundary):
            if target_kmer[4] == 'b':
                hashable_boundary = hashable_boundaries(target_boundary)
                matched_masses_b[hashable_boundary].append(target_kmer)
                kmer_index = kmer_index + 1

            if target_kmer[4] == 'y':
                hashable_boundary = hashable_boundaries(target_boundary)
                matched_masses_y[hashable_boundary].append(target_kmer)
                kmer_index = kmer_index + 1
            
        elif target_kmer[0] < target_boundary[0]:
            kmer_index = kmer_index + 1
            starting_point = starting_point + 1
        else:                                            #target_kmer > target_boundary[1]
            boundary_index = boundary_index + 1
            kmer_index = starting_point

    return matched_masses_b, matched_masses_y

# def modified_add_matched_to_matched_set(matched_masses_b_batch,matched_masses_b,kmer_set,batch_kmer_set,matched_masses_y_batch,matched_masses_y):
#     for k, v in matched_masses_b_batch.items():
#         matched_masses_b[k] += v 
#         for kmer in v:
#             kmer_set[kmer] += batch_kmer_set[kmer]
#     for k, v in matched_masses_y_batch.items():
#         matched_masses_y[k] += v 
#         for kmer in v:
#             kmer_set[kmer] += batch_kmer_set[kmer]

def modified_match_masses_per_protein(kv_prots,max_len,boundaries,kmer_set):
    extended_kv_prots = [(k, entry) for (k, v) in kv_prots for entry in v]
    kmers, kmer_set = modified_make_database_set(extended_kv_prots, max_len)
    # check_for_y_kmers(kmers)
    matched_masses_b, matched_masses_y = modified_merge(kmers, boundaries)
    # modified_add_matched_to_matched_set(matched_masses_b,kmer_set,kmers,matched_masses_y)

    return matched_masses_b, matched_masses_y, kmer_set

def modified_match_masses(boundaries: dict, db: Database, max_pep_len: int):
    # matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list) #Not sure this is needed
    max_boundary = max(boundaries.keys())
    estimated_max_len = ceil(boundaries[max_boundary][1] / 57.021464)
    max_len = min(estimated_max_len, max_pep_len)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    matched_masses_b, matched_masses_y, kmer_set = modified_match_masses_per_protein(kv_prots,max_len,boundaries,db)
    return (matched_masses_b, matched_masses_y, kmer_set)

def get_midpoints_of_boundries(spectra_boundaries):
    midpoints = []
    for spectra_boundry in spectra_boundaries:
        lower,upper = spectra_boundry
        midpoint = lower + ((upper-lower)/2)
        midpoints.append(midpoint)
    return midpoints        

def is_kmers_in_kmer_matches(kmers,kmer_matches):
    contained = False
    for kmer_match in kmer_matches:
        if kmer_match == kmers:
            contained = True
            break
    return contained

def get_speactras_for_mz_value(ion_charge, mz_value, ppm_tolerance):
    base_url = "http://hypedsearchservice.azurewebsites.net/api/proteinmatch?"
    url = base_url + "ion_charge=" + ion_charge
    url += "&weight=" + str(mz_value)
    url += "&ppm_tolerance=" + str(ppm_tolerance)
    request = requests.get(url = url)
    protein_matches = request.json()
    kmer_matches = []
    proteins_matched = []
    for protein_match in protein_matches:
        protein_match_list = list(protein_match.values())
        kmers = protein_match_list[1]
        protein = protein_match_list[0]
        proteins_matched.append(protein)
        if is_kmers_in_kmer_matches(kmers,kmer_matches):
            pass
        else:
            kmer_matches.append(kmers)
    return kmer_matches, proteins_matched

def match_masses_using_webservice(spectra_boundaries, ppm_tolerance):
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)
    midpoints = get_midpoints_of_boundries(spectra_boundaries)
    adjusted_ppm_tolerance = abs(ppm_tolerance / 1000)
    for index, midpoint in enumerate(midpoints):
        ion_charge = "B"
        boundry = str(spectra_boundaries[index][0]) + '-' + str(spectra_boundaries[index][1])
        matched_spectra,proteins_matched = get_speactras_for_mz_value(ion_charge, midpoint,adjusted_ppm_tolerance)
        matched_masses_b[boundry].append(matched_spectra)
        kmer_set = proteins_matched
    for index, midpoint in enumerate(midpoints):
        ion_charge = "Y"
        boundry = str(spectra_boundaries[index][0]) + '-' + str(spectra_boundaries[index][1])
        matched_spectra,proteins_matched = get_speactras_for_mz_value(ion_charge, midpoint,adjusted_ppm_tolerance)
        kmer_set.extend(proteins_matched)
        matched_masses_y[boundry].append(matched_spectra)
    return (matched_masses_b, matched_masses_y, kmer_set)
