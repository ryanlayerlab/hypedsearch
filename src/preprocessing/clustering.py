import collections
import operator
import os
from utils import ppm_to_da, to_percent
import gen_spectra
from constants import WATER_MASS, PROTON_MASS
from sqlite import database_file
import time
from constants import AMINO_ACIDS
from scoring.scoring import calc_bayes_score

def create_clusters_from_foos(foos):
    if len(foos) == 0 : return None
    clusters = []
    clusters.append(len(foos))
    clusters.append(foos[0].pid)
    max_len = 0
    max_hit = None
    for foo in foos:
        l = foo.end - foo.start + 1
        if l > max_len:
            max_len = l
            max_hit = foo
    clusters.append(max_hit.mz)
    clusters.append(max_hit.start)
    clusters.append(max_hit.end)
    clusters.append(max_hit.charge)
    clusters.append(foos)
    return clusters

def parse_foos(Hit, all_hits):
    hits = []
    #Matched masses data is of form (mass, start, end, ion_int, charge, protein_num)
    for A in all_hits:
        pid = int(A[2][5])
        start = int(A[2][1])
        end = int(A[2][2])
        mz = A[1]
        charge = A[2][4]

        hits.append( Hit(pid=pid, start=start, end=end, mz=mz, charge=charge) )
    return hits

def create_clusters(ion, b_hits, y_hits):
    all_clusters = []
    Foo = collections.namedtuple('Foo', 'pid start end mz charge')
    if ion == 'b':
        foos = parse_foos(Foo, b_hits)
        sorted_foos = sorted(foos, key=operator.attrgetter('pid', 'start', 'end'))
        last_pid = None
        last_start = None
        foos = []
        for foo in sorted_foos:
            if last_pid == foo.pid and last_start == foo.start:
                foos.append(foo)
            else:
                if foos != []:
                    clusters = create_clusters_from_foos(foos)
                    all_clusters.append(clusters)
                foos = [foo]
            last_pid = foo.pid
            last_start = foo.start
    else:
        foos = parse_foos(Foo, y_hits)
        sorted_foos = sorted(foos, key=operator.attrgetter('pid', 'end', 'start'))
        last_pid = None
        last_start = None
        foos = []
        for foo in sorted_foos:
            if last_pid == foo.pid and last_end == foo.end:
                foos.append(foo)
            else:
                if foos != []:
                    clusters = create_clusters_from_foos(foos)
                    all_clusters.append(clusters)
                foos = [foo]
            last_pid = foo.pid
            last_end = foo.end
    return all_clusters

def parse_indices(index_set):
    indices = []
    for index in index_set:
        string = str(index)
        A = string.rstrip().split(',')
        start = A[0]
        end = A[1]
        seq = A[2]
        mz = A[3]
        disallowed_characters = " ()\'"
        for character in disallowed_characters:
            start = start.replace(character, "")
            end = end.replace(character, "")
            seq = seq.replace(character, "")
            mz = mz.replace(character, "")
        target_tuple = (int(start), int(end), seq, float(mz))
        indices.append(target_tuple)
    return indices

def find_sequence(pid, start_ind, end_ind, proteins):
    protein = proteins[pid]
    prot_seq = protein[1]
    target = prot_seq[start_ind: end_ind]
    return target

def append_AA(next_AA, current_mass, ion, charge):
    
    raw_current_mass = gen_spectra.get_raw_mass(current_mass, ion, charge)
    new_raw = raw_current_mass + AMINO_ACIDS[next_AA]
    normalized_raw = gen_spectra.calc_combined_mass(new_raw, ion)
    return normalized_raw

def test_digest_match(protein_list, pid, digest, start, end, ion):
    parent_seq = protein_list[pid][1]
    if ion == 0:
        if parent_seq[start] == digest[0]: #if we correctly cut to the left. Means sequence starts with D
            return True
        target_position = start - 1
        if not target_position < 0:
            if parent_seq[target_position] == digest[1]: #if we correctly cut to the right. Means sequence starts after an R
                return True
    else:
        if parent_seq[end-1] == digest[1]: #if we correctly cut to the right. Means sequence ends in an R
            return True
        target_position = end
        if target_position < len(parent_seq):
            if parent_seq[target_position] == digest[0]: #if we correctly cut to the left. Means sequence ends to the left of a D
                return True
    
    return False

def find_extensions(conv_prec,current_mass,ion,charge,pid,protein_list,start,end,ppm_tolerance,seq,score):
    #goal is to get a list of all database extensions leading up to the precursor mass
    prot_seq = protein_list[pid][1]
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    extensions = []
    repeat = True
    current_score = score
    current_seq = seq
    
    if ion == 0:
        current_position = end
        while(repeat):
            if current_position < len(prot_seq):
                next_AA = prot_seq[current_position]
                if next_AA not in bad_chars:
                    current_mass = append_AA(next_AA, current_mass, ion, charge)
                    charge = 2
                    tol = ppm_to_da(current_mass, ppm_tolerance)
                    if current_mass < conv_prec + tol:
                        current_score += 1
                        current_seq = current_seq + next_AA
                        tup = (current_mass, start, current_position+1, 0, 2, pid, current_seq, current_score)
                        extensions.append(tup)
                        current_position += 1
                    else:
                        repeat = False
                else:
                    repeat = False
            else:
                repeat = False
    
    else:
        current_position = start-1
        while(repeat):
            if current_position >= 0:
                next_AA = prot_seq[current_position]
                if next_AA not in bad_chars:
                    current_mass = append_AA(next_AA, current_mass, ion, charge)
                    charge = 2
                    tol = ppm_to_da(current_mass, ppm_tolerance)
                    if current_mass < conv_prec + tol:
                        current_score += 1
                        current_seq = next_AA + current_seq
                        tup = (current_mass, current_position, end, 1, 2, pid, current_seq, current_score)
                        extensions.append(tup)
                        current_position = current_position - 1
                    else:
                        repeat = False
                else:
                    repeat = False
            else:
                repeat = False
    
    return extensions

def check_unique(merges):
    for m in merges:
        if merges.count(m) > 1:
            print(m,"has a count of", merges.count(m))
            return True
    return False                

def convert_components(component_arr, ion, score, seq):
    Foo = collections.namedtuple('Foo', 'pid start end mz charge')
    converted_components = []
    #(current_mass, start, end, 1, 2, pid)
    components_rev = list(reversed(component_arr))
    i = 0
    if ion == 0:
        prev_component = list(reversed(component_arr))[0]
        for component in list(reversed(component_arr)):
            if i != 0:
                while component.end < prev_component.end - 1:
                    new_component = Foo(pid = prev_component.pid, start = prev_component.start, end = prev_component.end - 1, mz = gen_spectra.max_mass(seq[:len(seq)-i], 'b', prev_component.charge), charge=prev_component.charge)
                    converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[:len(seq)-i], score-i))
                    prev_component = new_component
                    i +=1
                
            converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[:len(seq)-i], score-i))
            prev_component = component
            i += 1
    else:
        prev_component = component_arr[0]
        for component in component_arr:
            if i != 0:
                while component.start > prev_component.start +1:
                    new_component = Foo(pid = prev_component.pid, start = prev_component.start + 1, end = prev_component.end, mz = gen_spectra.max_mass(seq[i:], 'y', prev_component.charge), charge=prev_component.charge)
                    converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[i:], score-i))
                    prev_component = new_component
                    i +=1
            
            converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[i:], score-i))
            prev_component = component
            i += 1

    return converted_components

def old_score_clusters(ion, clusters, conv_prec, protein_list, prec_charge, ppm_tol):
    sorted_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge components seq')
    cluster_dict = dict()
    for i, A in enumerate(clusters):
        pid = A[1]
        mz = A[2]
        start = A[3]
        end = A[4]
        charge = A[5]
        seq = find_sequence(pid, start, end, protein_list)
        # digest_match = test_digest_match(protein_list,pid,digest,start,end,ion) #this code doesn't work
        score = A[0]
        components = convert_components(A[6], ion, score, seq)
        extensions = find_extensions(conv_prec,mz,ion,charge,pid,protein_list,start,end,ppm_tol,seq,score)
        target_cluster = sorted_cluster(score=score, pid=pid, start=start, end=end, mz=mz, charge=charge, components=components + extensions, seq=seq)
        converted_precursor = gen_spectra.convert_ion_to_precursor(mz, ion, charge, prec_charge)
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
        precursor = gen_spectra.convert_ion_to_precursor(b_mass, 0, b_charge, prec_charge)
    else: #y overlaps b
        precursor = gen_spectra.convert_ion_to_precursor(y_mass, 1, y_charge, prec_charge)
    return precursor

def calc_from_sequences(start, y_end, pid, max_len, prec_charge):
    db = database_file(max_len, False)
    entries = db.query_sequence(pid, start, y_end)
    if entries == []:
        return 0
    else:
        entry = entries.pop()
        precursor = gen_spectra.convert_ion_to_precursor(entry[0], entry[3], entry[4], prec_charge)
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
            combined_precursor = gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, b_charge, y_charge, precursor_charge)
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

# def filter_by_structure(natural_merges, hybrid_merges):
#     for merge in natural_merges:
#         pid = merge[3][0]
#         b_start, b_end = comb_seq[3][1], comb_seq[3][2]
#         y_start, y_end = comb_seq[4][1], comb_seq[4][2]
#         b_charge, y_charge = comb_seq[3][5], comb_seq[4][5]
#         b_score, y_score = comb_seq[3][3], comb_seq[3][4]
#         b_mass = comb_seq[3][4]
#         y_mass = comb_seq[4][4]
        
#         if (b_start )

# def combine_merges(pure_seqs, hybrid_seqs, target_num): #TODO
#     merged_top = []
#     pure_index, hybrid_index = 0, 0
#     while len(merged_top) < target_num and pure_index < len(pure_seqs) and hybrid_index < len(hybrid_seqs):
#         pure = pure_seqs[pure_index]
#         hybrid = hybrid_seqs[hybrid_index]
#         if pure[0] >= hybrid[0]: #We give ties to the non-hybrid sequences
#             merged_top.append(pure)
#             pure_index += 1
#         else:
#             merged_top.append(hybrid)
#             hybrid_index += 1
#     while len(merged_top) < target_num and pure_index < len(pure_seqs):
#         merged_top.append(pure_seqs[pure_index])
#         pure_index += 1
#     while len(merged_top) < target_num and hybrid_index < len(hybrid_seqs):
#         merged_top.append(hybrid_seqs[hybrid_index])
#         hybrid_index += 1

#     return merged_top

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
    #Given a cluster we want to find everything that it can pair with
    # It can pair with anything up to a certain mass 
    matches = []
    for key in indexed_clusters.keys():
        if key<=target_val: #if key is a valid key
            matches.append(key)
    return matches
    
# def index_by_precursor_mass(sorted_clusters, pc, ion):
#     indexed = dict()
#     for mz, charge in sorted_clusters.keys():
#         converted_precursor = gen_spectra.convert_ion_to_precursor(mz,ion,charge,pc) #TODO: C
#         if converted_precursor not in indexed.keys():
#             indexed[converted_precursor] = []
#         indexed[converted_precursor].append(mz)
#     indexed = collections.OrderedDict(sorted(indexed.items(),key=lambda t: t[0]))
#     return indexed
    
def get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, obs_prec, precursor_tol, prec_charge):
    merged_seqs = dict()
    #want clusters to be a dictionary where the keys are the input mz value this was matched to
    # ind_y = index_by_precursor_mass(y_sorted_clusters, prec_charge, 1) #Currently no functionality for overlap
    
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

def get_search_space(b_sorted_clusters, y_sorted_clusters, prec_charge): #This will eventually be reworked throughout the entire codebase but this is for proof of concept
    b_searches, y_searches = dict(), dict()
    for key in b_sorted_clusters.keys():
        for b in b_sorted_clusters[key]:
            for component in b.components:
                mass = component[0]
                charge = component[4]
                prec = gen_spectra.convert_ion_to_precursor(mass, 0, charge, prec_charge)
                if prec not in b_searches.keys():
                    b_searches[prec] = []
                b_searches[prec].append(component)
            
    for key in y_sorted_clusters.keys():       
        for y in y_sorted_clusters[key]:
            for component in y.components:
                mass = component[0]
                charge = component[4]
                prec = gen_spectra.convert_ion_to_precursor(mass, 1, charge, prec_charge)
                if prec not in y_searches.keys():
                    y_searches[prec] = []
                y_searches[prec].append(component)
            
    return b_searches, y_searches