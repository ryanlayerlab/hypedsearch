from utils import hashable_boundaries, predicted_len
from objects import Database
from collections import defaultdict
from typing import Iterable
from math import ceil

import gen_spectra
import array as arr

BATCH_SIZE = 300

def merge(
    mz_s: Iterable, 
    indices: Iterable, 
    kmers: Iterable, 
    boundaries: Iterable
    ) -> defaultdict:
    '''Perform a linear search of observed mz values that fit into mappings 
    to get a mapping from mz boundaries [*lower_bound*, *upper_bound*] to a set 
    of k-mers

    :param mz_s: mz values to look through
    :type mz_s: Iterable
    :param indices: index mappings from mz values to kmers. The steps to get 
        from an *m/z* value to a set of k-mers is as follows: if we have a m/z
        value at index *i*, we will get the values in range of *indices*[*i*-1] 
        to *indices*[*i*], call it j in J. Then, the k-mers we want are all kmers
        at *kmers*[*j*] for each j in J.
    :type indices: Iterable
    :param kmers: the k-mers associated with the mz_s in the range of indices
        as described in the *indices* param
    :type kmers: Iterable
    :param boundaries: lower upper bounds for a mass in a list like 
        [*lower_bound*, *upper_bound*]
    :type boundaries: Iterable

    :returns: mapping from a string (using src.utils.hashable_boundaries) 
        of <lower_bound>-<upper_bound> to a list of k-mers
    :rtype: defaultdict
    '''

    boudaries_index, spectra_index = 0, 0

    matched_masses = defaultdict(list)

    while boudaries_index < len(boundaries) and spectra_index < len(mz_s):

        # if mz_s[mz_i] is in the boundary, keep track of it increment mz_i
        if boundaries[boudaries_index][0] <= mz_s[spectra_index] <= boundaries[boudaries_index][1]:
            matched_masses[hashable_boundaries(boundaries[boudaries_index])] += kmers[indices[spectra_index - 1]:indices[spectra_index]]
            spectra_index += 1

        # if the upper boundary is less than the mz_i, increment b_i
        elif mz_s[spectra_index] > boundaries[boudaries_index][1]:
            boudaries_index += 1

        # if the lower boundary is greater than the mz_i, increment mz_i
        elif mz_s[spectra_index] < boundaries[boudaries_index][0]:
            spectra_index += 1

    return matched_masses

#todo: replace this with the web call
def make_database_set(
    proteins: list, 
    max_len: int
    ) -> (arr.array, arr.array, list, arr.array, arr.array, list, dict):
    '''Create parallel lists of (masses, index_maps, kmers) for the merge sort operation
    where index_maps map the massses to a range of positions in the kmers list 

    :param proteins: protein entries of shape (name, entry) where entry has a 
        *.sequence* attribute
    :type proteins: list
    :param max_len: max k-mer length
    :type max_len: int
    
    :returns: b ion masses created from the protein set, mapping from b ion masses 
        to the kmers associated (same length as b ion masses), kmers associated with 
        b ion masses, y ion masses created from the protein set, mapping from y ion masses 
        to the kmers associated (same length as y ion masses), kmers associated with 
        y ion masses, mapping of k-mers to source proteins
    :rtype: (array, array, list, array, array, list, dict)
    '''
    db_dict_b = defaultdict(set)
    db_dict_y = defaultdict(set)
    kmer_set = defaultdict(list)

    # function to add all masses of b+, b++, y+, y++ ions to the correct dict for sorting
    # and keep track of the source proteins
    def add_all(kmer, prot_name):
        for ion in 'by':
            for charge in [1, 2]:
                #hack
                #spec = pre_spec.spectrum
                #the_type = type(spec) #python = 'dict' #cpp = 'list'
                pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
                spec = pre_spec
                if isinstance(pre_spec,dict):
                    spec = pre_spec.get('spectrum')

                for i, mz in enumerate(spec):
                    kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
                    r_d = db_dict_b if ion == 'b' else db_dict_y
                    r_d[mz].add(kmer_to_add)
                    kmer_set[kmer_to_add].append(prot_name)

    plen = len(proteins)

    # todo: this is it
    # go through each protein and add all kmers to the correct dictionary for later sorting
    for i, (prot_name, prot_entry) in enumerate(proteins):
        
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        
        for j in range(1, max_len):
            kmer = prot_entry.sequence[:j]
            add_all(kmer, prot_name)

        for j in range(len(prot_entry.sequence) - max_len):
            kmer = prot_entry.sequence[j:j+max_len]
            add_all(kmer, prot_name)

        for j in range(len(prot_entry.sequence) - max_len, len(prot_entry.sequence)):
            kmer = prot_entry.sequence[j:]
            add_all(kmer, prot_name)

    print('\nSorting the set of protein masses...')
    
    # arrays take less space than lists
    db_list_b, index_list_b, kmer_list_b = arr.array('f'), arr.array('i'), []
    db_list_y, index_list_y, kmer_list_y = arr.array('f'), arr.array('i'), []

    # sort the keys to make sure that we are going through masses the correct way
    sorted_keys = sorted(db_dict_b.keys())
    for mz in sorted_keys:

        # get all kmers associated with this mass, append it the b ion list, and keep track 
        # of the kmers in the kmer list
        kmers = db_dict_b[mz]
        db_list_b.append(mz)
        offset = 0 if not len(index_list_b) else index_list_b[-1]
        index_list_b.append(len(kmers) + offset)
        kmer_list_b += kmers

    sorted_keys = sorted(db_dict_y.keys())
    for mz in sorted_keys:

        # get all kmers associated with this mass, append it the y ion list, and keep track 
        # of the kmers in the kmer list
        kmers = db_dict_y[mz]
        db_list_y.append(mz)
        offset = 0 if not len(index_list_y) else index_list_y[-1]
        index_list_y.append(len(kmers) + offset)
        kmer_list_y += kmers

    print('Done')

    return db_list_b, index_list_b, kmer_list_b, db_list_y, index_list_y, kmer_list_y, kmer_set

def match_masses(
    spectra_boundaries: list, 
    db: Database, 
    max_pep_len: int = 30
    ) -> (dict, dict, Database):

    matched_masses_b  = defaultdict(list)
    matched_masses_y= defaultdict(list)
    kmer_set = defaultdict(list)

    # estimate the max len
    estimated_max_len = ceil(spectra_boundaries[-1][1] / 57.021464)
    max_len = min(estimated_max_len, max_pep_len)

    #TODO: DB
    # calc the number of batches needed
    num_batches = ceil(len(db.proteins) / BATCH_SIZE)

    #TODO: DB
    # create batches of proteins in the form of (prot name, prot entry)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

    # go through each batch set and create the list representation, merge, and keep good prots
    for batch_num, batch_set in enumerate(batched_prots):

        print(f'On batch {batch_num + 1}/{num_batches}\n', end='')
        #todo start here
        #I think extended_batch_set is a list of protins (the database)
        extended_batch_set = [(k, entry) for (k, v) in batch_set for entry in v]

        # create our list representation
        batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = make_database_set(extended_batch_set, max_len)

        # find tha batch matched masses for both b and y ions
        matched_masses_b_batch = merge(batch_b_list, index_list_b, batch_kmer_b, spectra_boundaries)
        matched_masses_y_batch = merge(batch_y_list, index_list_y, batch_kmer_y, spectra_boundaries)

        # add these these hits to our function scoped set of masses and to the kmer set
        for k, v in matched_masses_b_batch.items():
            matched_masses_b[k] += v 
            
            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

        for k, v in matched_masses_y_batch.items():
            matched_masses_y[k] += v 

            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

    #update kmers in db to the kmer set
    db = db._replace(kmers=kmer_set)

    return (matched_masses_b, matched_masses_y, db)