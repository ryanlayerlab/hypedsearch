import collections
import operator
import os
from utils import ppm_to_da
from gen_spectra import get_precursor
from constants import WATER_MASS, PROTON_MASS
from scoring.scoring import calc_bayes_score

def write_cluster(cluster):
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
    O.append(max_hit.seq)
    O.append(max_hit.mz)
    O.append(max_hit.start)
    O.append(max_hit.end)
    O.append(max_hit.o_num)
    for hit in cluster:
        O.append( (hit.start, hit.end, hit.seq, hit.mz, hit.o_num) ) 
    return O

def parse_hits(Hit, all_hits):
    hits = []
    for A in all_hits:
        pid = int(A[2][1])
        start = int(A[2][3].split('-')[0])
        end = int(A[2][3].split('-')[1])
        seq = A[2][2]
        mz = A[1]
        occurance_num = A[3]

        hits.append( Hit(pid=pid, start=start, end=end, seq=seq, mz=mz, o_num = occurance_num) )
    return hits

def create_clusters(ion, b_hits, y_hits):
    clusters = []
    Hit = collections.namedtuple('Hit', 'pid start end seq mz o_num')
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

def Bayes_Score_clusters(ion, clusters, kmer_set):
    cluster = collections.namedtuple('cluster', 'prob score pid start end seq mz indices')
    if ion == 'b':
        b_cluster_array = []
        for A in clusters:
            score = A[0]
            pid = int(A[1])
            seq = A[2]
            mz = float(A[3])
            start = int(A[4])
            end = int(A[5])
            indices = A[6:]
            prob = calc_bayes_score(seq, mz, indices, kmer_set)
            target_cluster = cluster(prob=prob, score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices)

            b_cluster_array.append(target_cluster)

        b_sorted_clusters = sorted(b_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return b_sorted_clusters
    else:
        y_cluster_array = []
        for A in clusters:
            score = A[0]
            pid = int(A[1])
            seq = A[2]
            mz = float(A[3])
            start = int(A[4])
            end = int(A[5])
            indices = A[6:]
            prob = calc_bayes_score(seq, mz, indices, kmer_set)
            target_cluster = cluster(prob=prob, score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices)
            y_cluster_array.append(target_cluster)

        y_sorted_clusters = sorted(y_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return y_sorted_clusters
    
def Score_clusters(ion, clusters):
    cluster = collections.namedtuple('cluster', 'score pid start end seq mz o_num indices')
    if ion == 'b':
        b_cluster_array = []
        for A in clusters:
            score = A[0]
            pid = int(A[1])
            seq = A[2]
            mz = float(A[3])
            start = int(A[4])
            end = int(A[5])
            o_num = int(A[6])
            indices = A[7:]
            target_cluster = cluster(score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, o_num=o_num, indices=indices)

            b_cluster_array.append(target_cluster)

        b_sorted_clusters = sorted(b_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return b_sorted_clusters
    else:
        y_cluster_array = []
        for A in clusters:
            score = A[0]
            pid = int(A[1])
            seq = A[2]
            mz = float(A[3])
            start = int(A[4])
            end = int(A[5])
            o_num = int(A[6])
            indices = A[7:]
            target_cluster = cluster(score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, o_num=o_num, indices=indices)
            y_cluster_array.append(target_cluster)

        y_sorted_clusters = sorted(y_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return y_sorted_clusters

def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.seq, cluster.o_num)

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

    for pid in B:
        if pid not in Y:
            continue

        sorted_B = sorted(B[pid], key=operator.attrgetter('pid', 'start', 'end'))
        sorted_Y = sorted(Y[pid], key=operator.attrgetter('pid', 'start', 'end'))

        for b in sorted_B:
            y_i = bsearch(b.start, sorted_Y)

            if y_i >= len(sorted_Y): break

            y = sorted_Y[y_i]

            while y_i < len(sorted_Y) and y.start - b.end < 10:
                y = sorted_Y[y_i]
                merge_seqs.append((b.score + y.score, b.end - y.start, y.end-b.start,min_info(b), min_info(y)))
                y_i += 1
    return merge_seqs

def filter_by_precursor(mseqs, obs_prec, tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if b_seq != y_seq:
            new_seq = b_seq + y_seq
        else:
            new_seq = b_seq
        if not (get_precursor(new_seq, charge) > obs_prec + tol):
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
        y_start = comb_seq[4][1]
        b_end = comb_seq[3][2]
        if (y_start - b_end > 0) & (y_start - b_end < 10):
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
    for i, prot_name in enumerate(db.proteins):
        if i == target_prot:
            protein = db.proteins[prot_name]
            prot_seq = protein[0][1]
            to_add = prot_seq[target_index] if (target_index < len(prot_seq) and target_index > 0) else ''
            break
    
    return to_add

def filter_by_missing_mass(db, mseqs, obs_prec, tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        new_seq = overlap(comb_seq)
        dif = obs_prec + tol - get_precursor(new_seq, charge)
        if dif <= 1: #tol can vary but i'm not sure how much. Tol is .05 for spec 4 Other hacks are 2*tol
            filtered_seqs.append(comb_seq)
        else:
            next_b = modified_find_next_mass(comb_seq[3], 'b', db)
            b_seq = comb_seq[3][4]
            y_seq = comb_seq[4][4]
            b_dif = obs_prec + tol - get_precursor(b_seq + next_b + y_seq, charge)
            next_y = modified_find_next_mass(comb_seq[4], 'y', db)
            y_dif = obs_prec + tol - get_precursor(b_seq + next_y + y_seq, charge)
            if b_dif >= 0 or y_dif >= 0:
                filtered_seqs.append(comb_seq)
                
    return filtered_seqs

def combine_merges(pure_seqs, hybrid_seqs, target_num): #TODO
    merged_top = []
    pure_index, hybrid_index = 0,0
    if len(hybrid_seqs) == 0:
        return pure_seqs[:50]
    if len(pure_seqs) == 0:
        return hybrid_seqs[:50]
    while len(merged_top) < target_num:
        if len(pure_seqs) < pure_index:
            [merged_top.append(hybrid_seqs[x]) for x in hybrid_seqs[:min(target_num, len(hybrid_seqs))]]
            return merged_top
        if len(hybrid_seqs) < hybrid_index:
            [merged_top.append(pure_seqs[x]) for x in pure_seqs[:min(target_num, len(hybrid_seqs))]]
            return merged_top
        pure = pure_seqs[pure_index]
        hybrid = hybrid_seqs[hybrid_index]
        if pure[0] >= hybrid[0]: #We give ties to the non-hybrid sequences
            merged_top.append(pure)
            pure_index = pure_index + 1
        else:
            merged_top.append(hybrid)
            hybrid_index = hybrid_index + 1
    return merged_top

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

def grab_matches(b,indexed_clusters, target_val, ion):
    #Given a cluster we want to find everything that it can pair with
    # It can pair with anything up to a certain mass 
    current_index = 0
    matches = []
    for key in indexed_clusters.keys():
        if key<=target_val: #if key is a valid key
            for y in indexed_clusters[key]:
                if ion == 'b':
                    matches.append((b.score + y.score, b.end - y.start, y.end-b.start,min_info(b), min_info(y)))
                else:
                    matches.append((b.score + y.score, b.end - y.start, y.end-b.start,min_info(y), min_info(b)))
        else:
#             match, modified_seq = check_for_hybrid_overlap()
            break            
    return matches
    
def index_by_precursor_mass(sorted_clusters, pc):
    indexed = dict()
    for y in sorted_clusters:
        if get_precursor(y.seq, pc) not in indexed.keys():
            indexed[get_precursor(y.seq, pc)] = []
        indexed[get_precursor(y.seq, pc)].append(y)
    indexed = collections.OrderedDict(sorted(indexed.items(),key=lambda t: t[0]))
    return indexed
    
def get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, obs_prec, precursor_tol, charge):
    merged_seqs = []
    ind_b, ind_y = index_by_precursor_mass(b_sorted_clusters, charge),index_by_precursor_mass(y_sorted_clusters, charge)
    for i, cluster in enumerate(b_sorted_clusters[:10]):
        cluster_seq = cluster.seq
        cluster_mass = get_precursor(cluster_seq, charge)
        tol = ppm_to_da(obs_prec, precursor_tol)
        if not (cluster_mass > obs_prec + tol):
            diff = obs_prec + tol - cluster_mass + (charge * PROTON_MASS) + WATER_MASS
            merges = grab_matches(cluster,ind_y, diff, 'b')
            [merged_seqs.append(x) for x in merges]
    for i, cluster in enumerate(y_sorted_clusters[:10]):
        cluster_seq = cluster.seq
        cluster_mass = get_precursor(cluster_seq, charge)
        tol = ppm_to_da(obs_prec, precursor_tol)
        if not (cluster_mass > obs_prec + tol):
            diff = obs_prec + tol - cluster_mass + (charge * PROTON_MASS) + WATER_MASS
#             print(get_precursor(cluster_seq + 'DL', charge), obs_prec + tol)
            merges = grab_matches(cluster,ind_b, diff, 'y')
            [merged_seqs.append(x) for x in merges]

    merged_seqs.sort(key=lambda a: a[0], reverse=True)
    return merged_seqs

