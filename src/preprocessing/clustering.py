import collections
import operator

from itertools import groupby
import os
from lookups.utils import ppm_to_da, to_percent
import computational_pipeline.gen_spectra
from lookups.constants import WATER_MASS, PROTON_MASS
from preprocessing.sqlite_database import Sqllite_Database
import time
from lookups.constants import AMINO_ACIDS
from scoring.scoring import calc_bayes_score
from lookups.objects import ClusterItem, Cluster

def get_cluster_id(key):
    (protein_id, start_position) = key
    cluster_id = str(protein_id) + ':' + str(start_position)
    return cluster_id

def get_peptide(kmer, sqllite_database):
    protein_id = kmer.protein_id
    location_start = kmer.location_start
    location_end = kmer.location_end
    protein = sqllite_database.get_protein(protein_id)
    full_peptide = protein[2]
    peptide = full_peptide[location_start: location_end]
    return peptide

def get_full_peptide(protein_id, sqllite_database):
    protein = sqllite_database.get_protein(protein_id)
    full_peptide = protein[2]
    return full_peptide

def get_score(cluster_peptide,cluster_items):
    covered_positions = set()
    for cluster_item in cluster_items:
        kmer = cluster_item.kmer
        covered_positions.update(range(kmer.location_start, kmer.location_end + 1))
    protein_length = len(cluster_peptide)
    percentage_coverage = (len(covered_positions) / protein_length) * 100
    return percentage_coverage

def create_b_clusters(b_kmers,sqllite_database):
    all_clusters = []
    sorted_b_kmers = sorted(b_kmers, key=operator.attrgetter('protein_id', 'location_start'))
    for key, kmers in groupby(sorted_b_kmers, key=operator.attrgetter('protein_id', 'location_start')):
        cluster_items = []
        cluster_id = get_cluster_id(key)
        for kmer in kmers:
            peptide = get_peptide(kmer, sqllite_database)
            cluster_item = ClusterItem(key=cluster_id,kmer=kmer,peptide=peptide)
            cluster_items.append(cluster_item)
        (protein_id, start_position) = key
        full_peptide = get_full_peptide(protein_id, sqllite_database)
        score = get_score(full_peptide,cluster_items)
        cluster = Cluster(protein_id=protein_id,peptide=full_peptide,score=score,cluster_items=cluster_items)
        all_clusters.append(cluster)
    return all_clusters

def create_y_clusters(y_kmers,sqllite_database):
    all_clusters = []
    sorted_y_kmers = sorted(y_kmers, key=operator.attrgetter('protein_id', 'location_end'))
    for key, kmers in groupby(sorted_y_kmers, key=operator.attrgetter('protein_id', 'location_end')):
        cluster_items = []
        cluster_id = get_cluster_id(key)
        for kmer in kmers:
            peptide = get_peptide(kmer, sqllite_database)
            cluster_item = ClusterItem(key=cluster_id,kmer=kmer,peptide=peptide)
            cluster_items.append(cluster_item)
        (protein_id, start_position) = key
        full_peptide = get_full_peptide(protein_id, sqllite_database)
        score = get_score(full_peptide,cluster_items)
        cluster = Cluster(protein_id=protein_id,peptide=full_peptide,score=score,cluster_items=cluster_items)
        all_clusters.append(cluster)
    return all_clusters

def append_AA(next_AA, current_mass, ion, charge):
    raw_current_mass = computational_pipeline.gen_spectra.get_raw_mass(current_mass, ion, charge)
    new_raw = raw_current_mass + AMINO_ACIDS[next_AA]
    normalized_raw = computational_pipeline.gen_spectra.calc_combined_mass(new_raw, ion)
    return normalized_raw

def find_extensions(percursor,current_mass,ion,charge,pid,sqllite_database,start,end,ppm_tolerance,seq,score):
    protein = sqllite_database.get_protein(pid)
    prot_seq = protein[2]
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    extensions = []
    repeat = True
    current_seq = seq
    
    if ion == 0:
        current_mass, charge, repeat, current_seq = find_b_ion_extensions(percursor, current_mass, ion, charge, pid, start, end, ppm_tolerance, score, prot_seq, bad_chars, extensions, repeat, current_seq)
    else:
        find_y_ion_extension(percursor, current_mass, ion, charge, pid, start, end, ppm_tolerance, score, prot_seq, bad_chars, extensions, repeat, current_seq)
    
    return extensions

def find_y_ion_extension(conv_prec, current_mass, ion, charge, pid, start, end, ppm_tolerance, score, prot_seq, bad_chars, extensions, repeat, current_seq):
    current_position = start-1
    while(repeat):
        if current_position >= 0:
            next_AA = prot_seq[current_position]
            if next_AA not in bad_chars:
                current_mass = append_AA(next_AA, current_mass, ion, charge)
                charge = 2
                tol = ppm_to_da(current_mass, ppm_tolerance)
                if current_mass < conv_prec + tol:
                    current_seq = next_AA + current_seq
                    tup = (current_mass, current_position, end, 1, 2, pid, current_seq, score)
                    extensions.append(tup)
                    current_position = current_position - 1
                else:
                    repeat = False
            else:
                repeat = False
        else:
            repeat = False

def find_b_ion_extensions(conv_prec, current_mass, ion, charge, pid, start, end, ppm_tolerance, score, prot_seq, bad_chars, extensions, repeat, current_seq):
    current_position = end
    while(repeat):
        if current_position < len(prot_seq):
            next_AA = prot_seq[current_position]
            if next_AA not in bad_chars:
                current_mass = append_AA(next_AA, current_mass, ion, charge)
                charge = 2
                tol = ppm_to_da(current_mass, ppm_tolerance)
                if current_mass < conv_prec + tol:
                    current_seq = current_seq + next_AA
                    tup = (current_mass, start, current_position+1, 0, 2, pid, current_seq, score)
                    extensions.append(tup)
                    current_position += 1
                else:
                    repeat = False
            else:
                repeat = False
        else:
            repeat = False
    return current_mass,charge,repeat,current_seq

def convert_components(component_arr, ion, score, seq):
    Foo = collections.namedtuple('Foo', 'pid start end mz charge')
    converted_components = []
    components_rev = list(reversed(component_arr))
    i = 0
    if ion == 0:
        prev_component = list(reversed(component_arr))[0]
        for component in list(reversed(component_arr)):
            if component.end == prev_component.end:
                converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[:len(seq)-i], score-i))
                prev_component = component
            while component.end < prev_component.end:
                i +=1
                target_mz = computational_pipeline.gen_spectra.get_max_mass(seq[:len(seq)-i], 'b', prev_component.charge)
                new_component = Foo(pid = prev_component.pid, start = prev_component.start, end = prev_component.end - 1, mz = target_mz, charge=prev_component.charge)
                converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[:len(seq)-i], score-i))
                prev_component = new_component
    else:
        prev_component = component_arr[0]
        for component in component_arr:
            if component.start == prev_component.start:
                converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[i:], score-i))
                prev_component = component
            while component.start > prev_component.start +1:
                i +=1
                target_mz = computational_pipeline.gen_spectra.get_max_mass(seq[i:], 'y', prev_component.charge)
                new_component = Foo(pid = prev_component.pid, start = prev_component.start + 1, end = prev_component.end, mz = target_mz, charge=prev_component.charge)
                converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[i:], score-i))
                prev_component = new_component

    return converted_components

def old_score_clusters(ion, clusters, conv_prec, sqllite_database, precursor_charge, ppm_tolerance):
    sorted_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge components seq')
    cluster_dict = dict()
    for i, A in enumerate(clusters):
        pid = A[1]
        mz = A[2]
        start = A[3]
        end = A[4]
        charge = A[5]
        seq = find_sequence(pid, start, end, sqllite_database)
        score = A[0]
        components = convert_components(A[6], ion, score, seq)
        extensions = find_extensions(conv_prec,mz,ion,charge,pid,sqllite_database,start,end,ppm_tolerance,seq,score)
        target_cluster = sorted_cluster(score=score, pid=pid, start=start, end=end, mz=mz, charge=charge, components=components + extensions, seq=seq)
        converted_precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(mz, ion, charge, precursor_charge)
        if converted_precursor not in cluster_dict.keys():
            cluster_dict[converted_precursor] = []
        cluster_dict[converted_precursor].append(target_cluster)
        
    return cluster_dict
    
def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.mz, cluster.charge, cluster.components, cluster.seq)

def bsearch(key, Y):
        lo = -1
        hi = len(Y)
        mid = -1
        while (hi - lo > 1):
            mid = int((hi+lo) / 2)
            if Y[mid].start < key:
                lo = mid
            else:
                hi = mid
        return hi

def ysearch(key, B):
        lo = -1
        hi = len(B)
        mid = -1
        while (hi - lo > 1):
            mid = int((hi+lo) / 2)
            if B[mid].start < key:
                lo = mid
            else:
                hi = mid
        return hi

def Ryan_merge(b_sorted_clusters, y_sorted_clusters):
    merge_seqs = list()

    B = {}
    for mz in b_sorted_clusters.keys():
        for c in b_sorted_clusters[mz]:
            if c.pid not in B:
                B[c.pid] = []
            B[c.pid].append(c)

    Y = {}
    for mz in y_sorted_clusters.keys():
        for c in y_sorted_clusters[mz]:
            if c.pid not in Y:
                Y[c.pid] = []
            Y[c.pid].append(c)

    for i, pid in enumerate(B):
        if pid not in Y:
            continue

        sorted_B = sorted(B[pid], key=operator.attrgetter('pid', 'start', 'end'))
        sorted_Y = sorted(Y[pid], key=operator.attrgetter('pid', 'start', 'end'))
        
        for j, b in enumerate(sorted_B):
            y_i = bsearch(b.start, sorted_Y)

            if y_i >= len(sorted_Y): break

            y = sorted_Y[y_i]

            while y_i < len(sorted_Y) and y.start - b.end < 10:
                y = sorted_Y[y_i]
                merge_seqs.append((b.score + y.score, b, y))
                y_i += 1
                    
        for j, y in enumerate(sorted_Y):
            b_i = ysearch(y.start, sorted_B) 

            if b_i >= len(sorted_B): break

            b = sorted_B[b_i]

            while b_i < len(sorted_B) and y.start - b.end < 10:
                b = sorted_B[b_i]
                merge_seqs.append((b.score + y.score, b, y))
                b_i += 1
    return merge_seqs

def calc_from_total_overlap(side, b_mass, b_charge, y_mass, y_charge, prec_charge):
    if side: #b overlaps y
        precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(b_mass, 0, b_charge, prec_charge)
    else: #y overlaps b
        precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(y_mass, 1, y_charge, prec_charge)
    return precursor

def calc_from_sequences(start, y_end, pid, max_len, prec_charge):
    sqllite_database = Sqlite_Database(max_len, False)
    entries = sqllite_database.query_sequence(pid, start, y_end)
    if entries == []:
        return 0
    else:
        entry = entries.pop()
        precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(entry[0], entry[3], entry[4], prec_charge)
        return precursor
    
def total_overlap(b_pid, y_pid, b_start, y_start, b_end, y_end):
    if b_pid == y_pid:
        if b_start == y_start:
            if b_end <= y_end:
                return True, False
        if b_end == y_end:
            if b_start <= y_start:
                return True, True
    return False, False

def filter_by_precursor(mseqs, obs_prec, tol, precursor_charge, max_len):
    filtered_seqs = []
    for comb_seq in mseqs:
        b_pid, y_pid = comb_seq[1][0], comb_seq[2][0]
        b_start, b_end = comb_seq[1][1], comb_seq[1][2]
        y_start, y_end = comb_seq[2][1], comb_seq[2][2]
        b_charge, y_charge = comb_seq[1][5], comb_seq[2][5]
        b_mass, y_mass = comb_seq[1][4], comb_seq[2][4]
        full, side = total_overlap(b_pid, y_pid, b_start, y_start, b_end, y_end)
        if full:
            combined_precursor = calc_from_total_overlap(side, b_mass, b_charge, y_mass, y_charge, precursor_charge)
        elif b_start <= y_start and b_end <= y_end and y_start < b_end:
            combined_precursor = calc_from_sequences(b_start, y_end, b_pid, max_len, precursor_charge)
        else:
            combined_precursor = computational_pipeline.gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, b_charge, y_charge, precursor_charge)
        if not (combined_precursor > obs_prec + tol):
            filtered_seqs.append(comb_seq)
    return filtered_seqs

def get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start):
    seq = ''
    if y_start > b_end:
        return b_seq + y_seq
    else:
        for i in range(b_start, y_start):
            seq = seq + b_seq[i]
        return seq
    
def overlap(comb_seq):
    b_seq = comb_seq[3][4]
    y_seq = comb_seq[4][4]
    b_pid = comb_seq[3][0]
    y_pid = comb_seq[4][0]
    if b_pid == y_pid:
        if b_seq == y_seq:
            return b_seq
        else:
            y_start = comb_seq[4][1]
            b_end = comb_seq[3][2]
            if (y_start - b_end > 0) and (y_start - b_end < 10):
                b_start = comb_seq[3][1]
                return get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start)
            else:
                return b_seq + y_seq
    else:
        return b_seq + y_seq

def modified_find_next_mass(cluster, ion, db):
    if ion == 'b':
        target_index = cluster[2] + 1
    else:
        target_index = cluster[1]-1
    target_prot = cluster[0]
    protein = db.proteins[target_prot]
    prot_seq = protein[1]
    to_add = prot_seq[target_index] if (target_index < len(prot_seq) and target_index > 0) else ''
    return to_add

def check_for_hybrid_overlap(b_seq, y_seq, ion):
    match = True
    if ion == 'b':
        for i, char in enumerate(b_seq):
            if char == y_seq[0]:
                k = 0
                for j in range(i, len(b_seq) + 1):
                    if b_seq[j] != y_seq[k]:
                        match = False
                        break
        if match == True:
            print('Match was true for', b_seq)
            modified_seq = b_seq[:i]
    else:
        for i, char in enumerate(y_seq):
            if char == y_seq[0]:
                k = 0
                for j in range(i, len(b_seq) + 1):
                    if b_seq[j] != y_seq[k]:
                        match = False
                        break
        if match == True:
            print('Match was true for', b_seq)
            modified_seq = b_seq[:i]
    return match, modified_seq

def grab_y_matches(indexed_clusters,target_val):
    matches = []
    for key in indexed_clusters.keys():
        if key<=target_val: #if key is a valid key
            matches.append(key)
    return matches
        
def get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, obs_prec, precursor_tol, prec_charge):
    merged_seqs = dict()
    tol = ppm_to_da(obs_prec, precursor_tol)
    for conv_prec in b_sorted_clusters.keys():
        if not (conv_prec > obs_prec + tol):
            diff = obs_prec + tol - conv_prec + (prec_charge * PROTON_MASS) + WATER_MASS
            merges = grab_y_matches(y_sorted_clusters, diff)
            merged_seqs[conv_prec] = []
            [merged_seqs[conv_prec].append(x) for x in merges]
    return merged_seqs

def distribute_merges(merges, b_sorted_clusters, y_sorted_clusters):
    merged_clusters = []
    i = 0
    for key in merges.keys():
        for y_conv in merges[key]:
            for b in b_sorted_clusters[key]:
                for y in y_sorted_clusters[y_conv]:
                    if b.score + y.score > 4:
                        if i < 100000000:
                            merged_clusters.append((b.score + y.score, b, y))
                            i = i + 1
                        else:
                            return merged_clusters

    return merged_clusters

def get_search_space(sorted_clusters,precursor_charge):
    b_sorted_clusters = sorted_clusters.b_sorted_clusters
    y_sorted_clusters = sorted_clusters.y_sorted_clusters
    b_searches, y_searches = dict(), dict()
    for key in b_sorted_clusters.keys():
        for b in b_sorted_clusters[key]:
            for component in b.components:
                mass = component[0]
                charge = component[4]
                prec = computational_pipeline.gen_spectra.convert_ion_to_precursor(mass, 0, charge, precursor_charge)
                if prec not in b_searches.keys():
                    b_searches[prec] = []
                b_searches[prec].append(component)
            
    for key in y_sorted_clusters.keys():       
        for y in y_sorted_clusters[key]:
            for component in y.components:
                mass = component[0]
                charge = component[4]
                prec = computational_pipeline.gen_spectra.convert_ion_to_precursor(mass, 1, charge, precursor_charge)
                if prec not in y_searches.keys():
                    y_searches[prec] = []
                y_searches[prec].append(component)
    search_space = objects.SearchSpace(b_search_space=b_searches,y_search_space=y_searches)
    return search_space