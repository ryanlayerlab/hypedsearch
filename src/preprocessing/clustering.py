import collections
import operator
import os
from utils import ppm_to_da, to_percent
import gen_spectra
from constants import WATER_MASS, PROTON_MASS
from sqlite import database_file
import time
from constants import AMINO_ACIDS

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

# def Bayes_Score_clusters(ion, clusters, kmer_set):
#     cluster = collections.namedtuple('cluster', 'prob score pid start end seq mz indices')
#     if ion == 'b':
#         b_cluster_array = []
#         for A in clusters:
#             score = A[0]
#             pid = int(A[1])
#             seq = A[2]
#             mz = float(A[3])
#             start = int(A[4])
#             end = int(A[5])
#             indices = A[6:]
#             prob = calc_bayes_score(seq, mz, indices, kmer_set)
#             target_cluster = cluster(prob=prob, score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices)

#             b_cluster_array.append(target_cluster)

#         b_sorted_clusters = sorted(b_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
#         return b_sorted_clusters
#     else:
#         y_cluster_array = []
#         for A in clusters:
#             score = A[0]
#             pid = int(A[1])
#             seq = A[2]
#             mz = float(A[3])
#             start = int(A[4])
#             end = int(A[5])
#             indices = A[6:]
#             prob = calc_bayes_score(seq, mz, indices, kmer_set)
#             target_cluster = cluster(prob=prob, score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices)
#             y_cluster_array.append(target_cluster)

#         y_sorted_clusters = sorted(y_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
#         return y_sorted_clusters

def find_sequence(pid, start_ind, end_ind, proteins):
    protein = proteins[pid]
    prot_seq = protein[1]
    target = prot_seq[start_ind: end_ind]
    return target

# def find_extensions(pid, start, end, mass, ion, protein_list, charge, prec, prec_charge):
#     prot_seq = protein_list[pid][1]
#     current_mass = gen_spectra.convert_ion_to_precursor(mass,ion,charge,prec_charge)
#     extensions = []
#     bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
#     if ion == 0: #b ion
#         raw_b_mass = gen_spectra.get_raw_mass(mass, ion, charge)
#         current_position = end
#         while current_mass < prec:
#             if current_position > len(prot_seq)-1:
#                 break
#             next_AA = prot_seq[current_position]
#             if not any (x in bad_chars for x in next_AA):
#                 raw_b_mass = raw_b_mass + AMINO_ACIDS[next_AA]
#                 normalized_mass = gen_spectra.calc_combined_mass(raw_b_mass, 0)
#                 current_mass = gen_spectra.convert_ion_to_precursor(normalized_mass, 0, 2, prec_charge)
#                 tup = (normalized_mass, start, current_position+1, 0, 2, pid)
#                 extensions.append(tup)
#                 current_position = current_position + 1
#             else:
#                 break
#     else: #y ion
#         raw_y_mass = gen_spectra.get_raw_mass(mass, ion, charge)
#         current_position = start-1
#         while current_mass < prec:
#             if current_position < 0:
#                 break
#             next_AA = prot_seq[current_position]
#             if not any (x in bad_chars for x in next_AA):
#                 raw_y_mass = raw_y_mass + AMINO_ACIDS[next_AA]
#                 current_mass = gen_spectra.calc_combined_mass(raw_y_mass, 1)
#                 tup = (current_mass, current_position, end, 1, 2, pid)
#                 extensions.append(tup)
#                 current_position = current_position - 1
#             else:
#                 break
#     return(extensions)

def append_AA(next_AA, current_mass, ion, charge):
    
    raw_current_mass = gen_spectra.get_raw_mass(current_mass, ion, charge)
    new_raw = raw_current_mass + AMINO_ACIDS[next_AA]
    normalized_raw = gen_spectra.calc_combined_mass(new_raw, ion)
    return normalized_raw

def find_extensions(conv_prec,current_mass,ion,charge,pid,protein_list,start,end,ppm_tolerance):
    #goal is to get a list of all database extensions leading up to the precursor mass
    prot_seq = protein_list[pid][1]
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    extensions = []
    repeat = True

    #get next AA
    #Append AA and do calculation
    #Repeat as long as calculation is less than precursor

    
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
                        tup = (current_mass, start, current_position+1, 0, 2, pid)
                        extensions.append(tup)
                        current_position += 1
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
                        tup = (current_mass, current_position, end, 1, 2, pid)
                        extensions.append(tup)
                        current_position -= 1
                    else:
                        repeat = False
    
    return extensions

def check_unique(merges):
    for m in merges:
        if merges.count(m) > 1:
            print(m,"has a count of", merges.count(m))
            return True
    return False

def old_score_clusters(ion, clusters, conv_prec, protein_list, prec_charge, ppm_tol):
    sorted_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge extensions seq')
    cluster_dict = dict()
    for A in clusters:
        score = A[0]
        pid = int(A[1])
        mz = float(A[2])
        start = int(A[3])
        end = int(A[4])
        charge = A[5]
        seq = find_sequence(pid, start, end, protein_list)
        extensions = find_extensions(conv_prec,mz,ion,charge,pid,protein_list,start,end,ppm_tol)
        target_cluster = sorted_cluster(score=score, pid=pid, start=start, end=end, mz=mz, charge=charge, extensions=extensions, seq=seq)
        converted_precursor = gen_spectra.convert_ion_to_precursor(mz, ion, charge, prec_charge)
        if converted_precursor not in cluster_dict.keys():
            cluster_dict[converted_precursor] = []
        cluster_dict[converted_precursor].append(target_cluster)
        
    return cluster_dict
    
def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.mz, cluster.charge, cluster.extensions, cluster.seq)

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
    for key in merges.keys():
        for y_conv in merges[key]:
            for b in b_sorted_clusters[key]:
                for y in y_sorted_clusters[y_conv]:
                    merged_clusters.append((b.score + y.score, b, y))

    return merged_clusters