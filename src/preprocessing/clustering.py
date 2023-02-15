import collections
import operator
import os
from utils import ppm_to_da, to_percent
import gen_spectra
from constants import WATER_MASS, PROTON_MASS
from sqlite import database_file
import time
from constants import AMINO_ACIDS

def write_cluster(cluster):
    #returns cluster of the form (score, pid, mz, start, end)
    if len(cluster) == 0 : return None
    O = []
    O.append(len(cluster))
    O.append(cluster[0].pid)
    max_len = 0
    max_hit = None
    for hit in cluster:
        l = hit.end - hit.start + 1
        if l > max_len:
            max_len = l
            max_hit = hit
    O.append(max_hit.mz)
    O.append(max_hit.start)
    O.append(max_hit.end)
    O.append(max_hit.charge)

    return O

def parse_hits(Hit, all_hits):
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
    clusters = []
    Hit = collections.namedtuple('Hit', 'pid start end mz charge')
    if ion == 'b':
        hits = parse_hits(Hit, b_hits)
        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'start', 'end'))
        last_pid = None
        last_start = None
        cluster = []
        for hit in sorted_hits:
            if last_pid == hit.pid and last_start == hit.start:
                cluster.append(hit)
            else:
                if cluster != []:   
                    clusters.append(write_cluster(cluster))
                cluster = [hit]
            last_pid = hit.pid
            last_start = hit.start
    else:
        hits = parse_hits(Hit, y_hits)
        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'end', 'start'))
        last_pid = None
        last_start = None
        cluster = []
        for hit in sorted_hits:
            if last_pid == hit.pid and last_end == hit.end:
                cluster.append(hit)
            else:
                if cluster != []:
                    clusters.append(write_cluster(cluster))
                cluster = [hit]
            last_pid = hit.pid
            last_end = hit.end
    return clusters


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

def find_extensions(pid, start, end, mass, ion, protein_list, charge, prec):
    prot_seq = protein_list[pid][1]
    current_mass = mass
    extensions = []
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    if ion == 0: #b ion
        raw_b_mass = gen_spectra.get_raw_mass(mass, ion, charge)
        current_position = end
        while current_mass < prec:
            if current_position > len(prot_seq)-1:
                break
            next_AA = prot_seq[current_position]
            if not any (x in bad_chars for x in next_AA):
                raw_b_mass = raw_b_mass + AMINO_ACIDS[next_AA]
                current_mass = gen_spectra.calc_combined_mass(raw_b_mass, 0)
                tup = (current_mass, start, current_position, 0, 2, pid)
                extensions.append(tup)
                current_position = current_position + 1
            else:
                break
    else: #y ion
        raw_y_mass = gen_spectra.get_raw_mass(mass, ion, charge)
        current_position = start-1
        while current_mass < prec:
            if current_position < 0:
                break
            next_AA = prot_seq[current_position]
            if not any (x in bad_chars for x in next_AA):
                raw_y_mass = raw_y_mass + AMINO_ACIDS[next_AA]
                current_mass = gen_spectra.calc_combined_mass(raw_y_mass, 1)
                tup = (current_mass, current_position, end, 1, 2, pid)
                extensions.append(tup)
                current_position = current_position - 1
            else:
                break
    return(extensions)

def get_all_extensions(all_side, ion, protein_list, prec):
    all_tuples = []
    for hit in all_side:
        score = hit[0]
        pid = hit[1]
        mz = hit[2]
        start = int(hit[3])
        end = int(hit[4])
        charge = hit[5]
        
        extensions = find_extensions(pid, start, end, mz, ion, protein_list, charge, prec)
        all_tuples.append((pid, start, end, score, mz, charge, extensions))
    return all_tuples

def get_unique_merged(merged, dict_b, dict_y, protein_list, b_prec, y_prec):
    #Problem with uniqueness at the moment since combinations can still have multiple pieces in common
    #Need to do extensions on each unique piece that shows up and then substitute them into existing merges    
    #Make a mapping between each unique cluster and all of it's extended replacements
    unique_left, unique_right = set(), set()
    b_indices, y_indices = dict(), dict()
    mb_dict = dict()
    for m in merged:
        unique_left.add(m[1])        
        unique_right.add(m[2])
        
    for m in unique_left:
        left_score = m[3]
        left_mz = m[4]
        left_len = m[2]-m[1]
        left_key = (left_score, left_len, left_mz)
        all_left = dict_b[left_key]
        extended_unique_left = get_all_extensions(all_left, 0, protein_list, b_prec)
        if m not in b_indices.keys():
            b_indices[m] = []
        [b_indices[m].append(x) for x in extended_unique_left]
        if m not in mb_dict.keys():
            mb_dict[m] = []
        mb_dict[m].append(m[2]) #m[2] is made from unique_left where we actually want the corresponding right. Does this live in any of our dicts?
    
    for m in unique_right:
        right_score = m[3]
        right_mz = m[4]
        right_len = m[2] - m[1]
        right_key = (right_score, right_len, right_mz)
        all_right = dict_y[right_key]
        extended_unique_right = get_all_extensions(all_right, 1, protein_list, y_prec)
        if m not in y_indices.keys():
            y_indices[m] = []
        [y_indices[m].append(x) for x in extended_unique_right]
        
    
    #now want a dictionary mapping each index to the unique place where they occur

    return b_indices, y_indices, mb_dict

def check_unique(merges):
    for m in merges:
        if merges.count(m) > 1:
            print(m,"has a count of", merges.count(m))
            return True
    return False

def extract_natives(large_merges):
    #Goal is to take in the full expanded list of merges and sort them into naturals or hybrids
    natural_merges, hybrid_merges = [], []
    for merge in large_merges:
        left_part = merge[3]
        right_part = merge[4]
        left_pid, right_pid = left_part[5], right_part[5]
        if left_pid == right_pid: #first requirement for native
            b_end = left_part[1]
            y_start = right_part[0]
            if y_start - b_end > 0 and y_start - b_end <= 10: #second requirement for native
                natural_merges.append(merge)
                continue
        hybrid_merges.append(merge)
        
    return natural_merges, hybrid_merges
                

def expand_clusters(merged, dict_b, dict_y, protein_list, b_prec, y_prec):
    merges = []    
    
    #key_tuple = (score, size, mz)
    for i,merge in enumerate(merged):
        left_part = merge[1]
        score = left_part[3]
        mz = left_part[2]
        size = left_part[1]-left_part[0]
        all_extended_left = dict_b[(score,size,mz)]
        
        right_part = merge[2]
        score = right_part[3]
        mz = right_part[2]
        size = right_part[1]-right_part[0]
        all_extended_right = dict_y[(score,size,mz)]        

        # print("all_extended_left keys are not unique:", check_keys_unique(all_extended_left.keys()))
        # print("all_extended_right keys are not unique:", check_keys_unique(all_extended_right.keys()))
        
        for left in all_extended_left: # These lines are not unique
            for right in all_extended_right:
                #left_score + right_score, b.end - y.start, y.end-b.start
                #score, pid, mz, start, end, charge
                lextensions = find_extensions(left[1], left[3], left[4], left[2], 0, protein_list, left[5], b_prec)
                rextensions = find_extensions(right[1], right[3], right[4], right[2], 1, protein_list, right[5], y_prec)
                nleft = (left[0], left[1], left[2], left[3], left[4], left[5], lextensions)
                nright = (right[0], right[1], right[2], right[3], right[4], right[5], rextensions)
                merges.append((left[3] + right[3], left[2]-right[1], right[2]-left[1], nleft, nright))
        
        # if check_unique(merges):
        #     print("Entry", i, "is not unique")
        # else:
        #     print("Entry", i, "is unique")
            
        
    return merges
    
def get_unique_clusters(clusters):
    # Want to determine uniqueness by the size and the mass matched to
    # Load them into dictionary
    cluster = collections.namedtuple('cluster', 'score pid start end mz charge')
    cluster_mapping = dict()
    unique_clusters = []
    for A in clusters:
        score = A[0]
        pid = int(A[1])
        mz = float(A[2])
        charge = A[5]
            
        start = int(A[3])
        end = int(A[4])

        size = end-start
        key_tuple = (score, size, mz)
        if key_tuple not in cluster_mapping.keys():
            cluster_mapping[key_tuple] = []
        cluster_mapping[key_tuple].append(A)
        
    for key in cluster_mapping:
        first_cluster = cluster_mapping[key][0]
        score = first_cluster[0]
        pid = first_cluster[1]
        mz = first_cluster[2]
        start = first_cluster[3]
        end = first_cluster[4]
        charge = first_cluster[5]
        
        target_cluster = cluster(score=score, pid=pid, start=start, end=end, mz=mz, charge=charge)

        unique_clusters.append(target_cluster)
        
    sorted_clusters = sorted(unique_clusters, key=operator.attrgetter('score', 'pid'), reverse = True)
        
    return cluster_mapping, sorted_clusters
    
def old_score_clusters(ion, clusters, conv_prec, protein_list):
    cluster = collections.namedtuple('cluster', 'score pid start end mz charge extensions')
    cluster_array = []
    # with open('clusters.txt', 'a') as b:
    for i, A in enumerate(clusters):
        score = A[0]
        pid = int(A[1])
        mz = float(A[2])
        start = int(A[3])
        end = int(A[4])
        charge = A[5]
        # query_time = time.time()
        extensions = find_extensions(pid, start, end, mz, ion, protein_list, charge, conv_prec)
        # extensions = []
        target_cluster = cluster(score=score, pid=pid, start=start, end=end, mz=mz, charge=charge, extensions=extensions)
        # print("Time:", time.time() - query_time, "End-Start", end-start)

        cluster_array.append(target_cluster)
        
        # b.write(str(score) + '\t' + str(pid) + '\t' + str(start) + '\t' + str(end) + '\t' + str(mz) + '\t' + str(charge) + '\t' + str(extensions) + '\n')

    sorted_clusters = sorted(cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)

    return sorted_clusters

def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.mz, cluster.charge, cluster.extensions)

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
    merge_seqs = []

    B = {}
    for c in b_sorted_clusters:
        if c.pid not in B:
            B[c.pid] = []
        B[c.pid].append(c)

    Y = {}
    for c in y_sorted_clusters:
        if c.pid not in Y:
            Y[c.pid] = []
        Y[c.pid].append(c)

    for i, pid in enumerate(B):
        if pid not in Y:
            continue

        sorted_B = sorted(B[pid], key=operator.attrgetter('pid', 'start', 'end'))
        sorted_Y = sorted(Y[pid], key=operator.attrgetter('pid', 'start', 'end'))
        
        for j, b in sorted_B:
            y_i = bsearch(b.start, sorted_Y)

            if y_i >= len(sorted_Y): break

            y = sorted_Y[y_i]

            while y_i < len(sorted_Y) and y.start - b.end < 10:
                y = sorted_Y[y_i]
                merge_seqs.append((b.score + y.score, min_info(b), min_info(y)))
                y_i += 1
                print(i,j)
                b_unique = check_unique(merge_seqs)
                    
        for j, y in sorted_Y:
            b_i = ysearch(y.start, sorted_B) 

            if b_i >= len(sorted_B): break

            b = sorted_B[b_i]

            while b_i < len(sorted_B) and y.start - b.end < 10:
                b = sorted_B[b_i]
                merge_seqs.append((b.score + y.score, min_info(b), min_info(y)))
                b_i += 1
                if b_unique != check_unique(merge_seqs):
                    print("y added something unique at", i, j)
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
        # checking cases
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

def grab_matches(b,indexed_clusters, target_val):
    #Given a cluster we want to find everything that it can pair with
    # It can pair with anything up to a certain mass 
    matches = []
    for key in indexed_clusters.keys():
        if key<=target_val: #if key is a valid key
            for y in indexed_clusters[key]:
                    matches.append((b.score + y.score, min_info(b), min_info(y)))
    return matches
    
def index_by_precursor_mass(sorted_clusters, pc, ion):
    indexed = dict()
    for y in sorted_clusters:
        converted_precursor = gen_spectra.convert_ion_to_precursor(y.mz,ion,y.charge,pc)
        if converted_precursor not in indexed.keys():
            indexed[converted_precursor] = []
        indexed[converted_precursor].append(y)
    indexed = collections.OrderedDict(sorted(indexed.items(),key=lambda t: t[0]))
    return indexed
    
def get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, obs_prec, precursor_tol, prec_charge):
    merged_seqs = []
    ind_b, ind_y = index_by_precursor_mass(b_sorted_clusters, prec_charge, 1), index_by_precursor_mass(y_sorted_clusters, prec_charge, 1) #Currently no functionality for overlap
    for cluster in b_sorted_clusters[:10]:
        cluster_mass = gen_spectra.convert_ion_to_precursor(cluster.mz, 0, cluster.charge, prec_charge)
        tol = ppm_to_da(obs_prec, precursor_tol)
        if not (cluster_mass > obs_prec + tol):
            diff = obs_prec + tol - cluster_mass + (prec_charge * PROTON_MASS) + WATER_MASS
            merges = grab_matches(cluster,ind_y, diff)
            [merged_seqs.append(x) for x in merges]
        # The commented code below doesn't need to be here if we are finding uniqueness from every b cluster
    for cluster in y_sorted_clusters[:10]: 
        cluster_mass = gen_spectra.convert_ion_to_precursor(cluster.mz, 1, cluster.charge, prec_charge)
        tol = ppm_to_da(obs_prec, precursor_tol)
        if not (cluster_mass > obs_prec + tol):
            diff = obs_prec + tol - cluster_mass + (prec_charge * PROTON_MASS) + WATER_MASS
#             print(get_precursor(cluster_seq + 'DL', charge), obs_prec + tol)
            merges = grab_matches(cluster,ind_b, diff)
            [merged_seqs.append(x) for x in merges]

    merged_seqs.sort(key=lambda a: a[0], reverse=True)
    return merged_seqs

